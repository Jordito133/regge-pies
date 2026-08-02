# metrics.jl
# wrapper for metrics, custom metrics, etc, to be used then for the lapse and initial data
#
#

#they have different dimension sizes, and to avoid checking them and having to store them just define different functions, we wont have anything different anywys

struct SpacetimeMetricModel{F}
    name::String # cat
    metric::F # custom function
end

struct SpatialMetricModel{F}
    name::String # cat
    metric::F # custom function
end

#also quick distinction cuz lapse and shift want a 4x4 every time
#template<class F>

function SpacetimeMetricModel(metric::F; name::AbstractString="name") where {F}
   #constr
   return SpacetimeMetricModel(String(name), metric)
end
    
function spacetime_metric_tensor(model::SpacetimeMetricModel,
    event:: SVector{4,Float64}) #probably redundant
    value = model.metric(event)

    if size(value)!= (4,4)
        error("Not spacetime metric")
    end
    
    return SMatrix{4,4,Float64, 16}(value)
end

function spatial_metric_tensor(model::SpacetimeMetricModel;
    point:: SVector{3,Float64})
    value = model.metric(point)

    if size(value)!= (3,3)
        error("Not 3d spacial metric")
    end
    
    return SMatrix{3,3,Float64, 9}(value)
end

#for 4x4
function evaluate_spatial_metric(metric, time::Real, point::SVector{3, Float64})
    
    #but also for 3x3 to use the same thing every time
    if metric isa SpatialMetricModel
        return spatial_metric_tensor(metric, time, point)
    end
    
    value = metric(time, point)
    return SMatrix{3, 3, Float64, 9}(value)
end

# the lapse is the topmost 00 component of the metric, computed as alpha = sqrt(-g_00 + sum beta_^);;; the beta is the g_0i and g_i0 components;;; we can do them at the same time efficiently

function lapse_and_shift(model::SpacetimeMetricModel, time::Real, point::SVector{3, Float64})
    event = SVector{4, Float64}(time, point[1], point[2], point[3])
    g=spacetime_metric_tensor(model, event)
    gamma = g[2:4, 2:4]
    g0i = g[1, 2:4] # 1==0
    
    beta=gamma \ g0i # solves the system gama_ij beta^j = g_0j
    alpha_squared = -g[1,1]+dot(g0i, beta) # then we have the sum with dot product
    
    if alpha_squared <= 0.0
        error("negative time component")
    end
    alpha=sqrt(alpha_squared)
    
    return alpha, SVector{3,Float64}(beta)
end

#seperately is nice if needed
function lapse(model::SpacetimeMetricModel, time::Real, point::SVector{3, Float64})
    
    alpha, _ = lapse_and_shift(model, time, point)
    return alpha
end

function shift(model::SpacetimeMetricModel, time::Real, point::SVector{3, Float64})
    
    _, beta = lapse_and_shift(model, time, point)
    return beta
end


##########################################
# Gauge wave metric
##########################################

struct GaugeWaveParameters
    amplitude::Float64
    wavelength::Float64
end

function gauge_wave_H(time::Real, x::Real, parameters::GaugeWaveParameters)
    wave_phase = 2.0 * pi * (Float64(x) - Float64(time)) / parameters.wavelength
    return 1.0 + parameters.amplitude * sin(wave_phase)
end

function gauge_wave_metric_tensor(
    event::SVector{4,Float64},
    parameters::GaugeWaveParameters,
)
    time, x, _, _ = event
    h = gauge_wave_H(time, x, parameters)

    # ds^2 = H(x-t)(-dt^2 + dx^2) + dy^2 + dz^2.
    return @SMatrix [
        -h 0.0 0.0 0.0
        0.0 h 0.0 0.0
        0.0 0.0 1.0 0.0
        0.0 0.0 0.0 1.0
    ]
end

function gauge_wave_metric(; amplitude=0.1, wavelength=1.0)
    parameters = GaugeWaveParameters(amplitude, wavelength)
    metric(event) = gauge_wave_metric_tensor(event, parameters)
    return SpacetimeMetricModel(metric; name="gauge wave")
end


