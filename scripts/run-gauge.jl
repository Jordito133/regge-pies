using ReggePIES
using Printf

function gauge_wave_x_edge_profile(template, edge_s)
    sums = zeros(Float64, ReggePIES.Nx)
    counts = zeros(Int, ReggePIES.Nx)

    for k in 0:(ReggePIES.Nz - 1), j in 0:(ReggePIES.Ny - 1), i in 0:(ReggePIES.Nx - 1)
        a = ReggePIES.vertex_id(i, j, k)
        b = ReggePIES.vertex_id(i + 1, j, k)
        eid = ReggePIES.edge_id(template, a, b)
        idx = i + 1
        sums[idx] += edge_s[eid]
        counts[idx] += 1
    end

    return [counts[i] == 0 ? NaN : sums[i] / counts[i] for i in eachindex(sums)]
end

template, coords, state = ReggePIES.build_gauge_wave_initial_state()
validation = ReggePIES.validate_gauge_wave_initial_data(template, coords, state.current.edge_s)
profile = gauge_wave_x_edge_profile(template, state.current.edge_s)
amp = 0.5 * (maximum(profile) - minimum(profile))

println("Gauge-wave initial metric validation")
println("  grid = ($(ReggePIES.Nx), $(ReggePIES.Ny), $(ReggePIES.Nz))")
println("  vertices = $(template.num_vert), edges = $(length(template.edges)), tetras = $(length(template.tetras))")
println("  checked_tetras = $(validation.checked_tetras)")
@printf("  max_gram_error = %.6e\n", validation.max_gram_error)
@printf("  max_metric_error = %.6e\n", validation.max_metric_error)
println("  worst_tetra = $(validation.worst_tetra)")
@printf("  x_edge_profile_min = %.8e\n", minimum(profile))
@printf("  x_edge_profile_max = %.8e\n", maximum(profile))
@printf("  x_edge_profile_amplitude = %.8e\n", amp)
println("  x_edge_profile = $(profile)")

