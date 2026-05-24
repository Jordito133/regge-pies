# geometry.jl
# Geometry formulas, related to combinatorics.jl
# Gram matrix

trianglekey_sorted(tri::NTuple{3,Int}) = trianglekey(tri...)

triangle_Q_from_sq(x::Real, y::Real, z::Real) = 2.0 * (x * y + x * z + y * z) - x^2 - y^2 - z^2

function edge_sq(edge_map::AbstractDict{NTuple{2,Int},<:Real}, i::Int, j::Int)
    edge = edgekey(i, j)
    haskey(edge_map, edge) || error("no lenght $edge.")
    return Float64(edge_map[edge])
end

# Gram matrix functions
function tetra_gram_matrix_from_sq(
    s01::Real, s02::Real, s03::Real,
    s12::Real, s13::Real, s23::Real,
)
    return [
        Float64(s01) 0.5 * (s01 + s02 - s12) 0.5 * (s01 + s03 - s13)
        0.5 * (s02 + s01 - s12) Float64(s02) 0.5 * (s02 + s03 - s23)
        0.5 * (s03 + s01 - s13) 0.5 * (s03 + s02 - s23) Float64(s03)
    ]
end

function tetra_gram_matrix_from_sq(tet::NTuple{4,Int}, edge_map::AbstractDict{NTuple{2,Int},<:Real})
    a, b, c, d = tet
    return tetra_gram_matrix_from_sq(
        edge_sq(edge_map, a, b),
        edge_sq(edge_map, a, c),
        edge_sq(edge_map, a, d),
        edge_sq(edge_map, b, c),
        edge_sq(edge_map, b, d),
        edge_sq(edge_map, c, d),
    )
end

# metric reconstruction

function reconstruct_metric_from_gram(gram::AbstractMatrix{<:Real}, edge_matrix::AbstractMatrix{<:Real})
    if size(gram, 1) != size(gram, 2)
        error("Square matrix required")
    end
    
    if size(edge_matrix, 2) == size(gram, 1)
        error("Same dimensions required")
    end
    
    inv_edge_matrix = inv(Matrix{Float64}(edge_matrix))
    return transpose(inv_edge_matrix) * Matrix{Float64}(gram) * inv_edge_matrix
end

