# state.jl
# Dynamics; states of the current and evolved slices


# a struct that includes the edge lengths
mutable struct SliceState
    slice_index::Int
    edge_s::Vector{Float64}
end

# temporary data that lives between the two slices
mutable struct SweepState
    #tester that tells us if a vertex has been evolved during this sweep
    evolved::Vector{Bool}
    #s_vv', initially empty
    vertical_s::Vector{Union{Nothing,Float64}}
    #s_v'u
    diagonal_s::Dict{Tuple{Int,Int},Float64}
end

mutable struct PIESState
    template::SliceTemplate #its crucial that the topology remains the same, and connections and graphing colors also remain the same through the evolution
    current::SliceState
    next::SliceState
    sweep::SweepState
end
# initial state
function initialize_state(template::SliceTemplate, current_edge_s::Vector{Float64}; slice_index::Int=0)
    #next slice is unkown
    next_edge_s = fill(NaN, length(template.edges))
    #also no known edges for the sweep
    sweep = SweepState(
        fill(false, template.num_vert),
        fill(nothing, template.num_vert),
        Dict{Tuple{Int,Int},Float64}(),
    )
    return PIESState(
        template,
        SliceState(slice_index, copy(current_edge_s)),
        SliceState(slice_index + 1, next_edge_s),
        sweep,
    )
end
#reset sweep data for next sweep
function start_new_sweep!(state::PIESState)
    fill!(state.sweep.evolved, false)
    fill!(state.sweep.vertical_s, nothing)
    empty!(state.sweep.diagonal_s)
    fill!(state.next.edge_s, NaN)
    state.next.slice_index = state.current.slice_index + 1
end
#promote next slice to current slice
# use start_new_sweep!()
function finish_sweep!(state::PIESState)
    state.current.edge_s .= state.next.edge_s
    state.current.slice_index = state.next.slice_index
    start_new_sweep!(state)
end


