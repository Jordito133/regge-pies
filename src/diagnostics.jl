# diagnostics.jl
#
# 

# validation  of gauge wave starting data
# G_{ij}=e_i^T \gamma e_j (in Xyz this is flat and just e_i^T e_j )
# so check if exact reconstruction leads to this

function validate_gauge_wave_initial_data(
    template::SliceTemplate,
    vertex_coords::Vector{SVector{3,Float64}},
    edge_s::AbstractVector{<:Real};
    max_tetras::Union{Nothing,Int}=nothing,
)
    edge_map = build_edge_sq_map(template, edge_s)
    Xvals, CX = build_gauge_wave_X_table()
    identity3 = Matrix{Float64}(I, 3, 3)

    gram_max_error = 0.0
    metric_max_error = 0.0
    worst_tetra = nothing

    last_tid = max_tetras === nothing ? length(template.tetras) : min(max_tetras, length(template.tetras))
    for tid in 1:last_tid
        tet = template.tetras[tid]
        base = vertex_coords[tet[1]]
        edge_matrix = hcat(
            wrapped_intrinsic_displacement(base, vertex_coords[tet[2]], Xvals, CX),
            wrapped_intrinsic_displacement(base, vertex_coords[tet[3]], Xvals, CX),
            wrapped_intrinsic_displacement(base, vertex_coords[tet[4]], Xvals, CX),
        )

        gram = tetra_gram_matrix_from_sq(tet, edge_map)
        gram_from_vectors = transpose(edge_matrix) * edge_matrix
        metric = reconstruct_metric_from_gram(gram, edge_matrix)

        local_gram_error = maximum(abs.(gram - gram_from_vectors))
        local_metric_error = maximum(abs.(metric - identity3))

        if max(local_gram_error, local_metric_error) > max(gram_max_error, metric_max_error)
            worst_tetra = tet
        end
        gram_max_error = max(gram_max_error, local_gram_error)
        metric_max_error = max(metric_max_error, local_metric_error)
    end

    return (
        checked_tetras = last_tid,
        max_gram_error = gram_max_error,
        max_metric_error = metric_max_error,
        worst_tetra = worst_tetra,
        passed = isfinite(gram_max_error) && isfinite(metric_max_error),
    )
end
