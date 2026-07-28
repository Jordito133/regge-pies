# geometry.jl
# Geometry formulas, related to combinatorics.jl
# Gram matrix

trianglekey_sorted(tri::NTuple{3,Int}) = trianglekey(tri...)

triangle_Q_from_sq(x::Real, y::Real, z::Real) = 2.0 * (x * y + x * z + y * z) - x^2 - y^2 - z^2

triangle_area_sq_from_sq(x::Real, y::Real, z::Real) =
    abs(triangle_Q_from_sq(x, y, z)) / 16.0

triangle_area_from_sq(x::Real, y::Real, z::Real) =
    sqrt(triangle_area_sq_from_sq(x, y, z))

# dA/dx = (y + z - x)/(16 A), using the positive area magnitude.
function d_triangle_area_d_sq(
    x::Real,
    y::Real,
    z::Real;
    null_tol::Real=1.0e-14,
)
    q = triangle_Q_from_sq(x, y, z)
    abs(q) > null_tol ||
        error("The triangle is null or degenerate.")
    
    return (y + z - x) / (16.0 * triangle_area_from_sq(x, y, z))
end

function build_edge_sq_map(template::SliceTemplate, edge_s::AbstractVector{<:Real})
    length(edge_s) == length(template.edges) || error("edge_s length does not match template edge count.")
    edge_map = Dict{NTuple{2,Int},Float64}()
    for (eid, edge) in enumerate(template.edges)
        edge_map[edge] = Float64(edge_s[eid])
    end
    return edge_map
end

function edge_sq(edge_map::AbstractDict{NTuple{2,Int},<:Real}, i::Int, j::Int)
    edge = edgekey(i, j)
    if !haskey(edge_map, edge)
        error("Missing squared length for edge $edge.")
    end
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

#4x4 spacetime gram matrix
function simplex_gram_matrix_from_sq(
    simplex::NTuple{5,Int},
    edge_map,
)
    reference = simplex[1]
    gram = Matrix{Float64}(undef, 4, 4)
    for a in 1:4, b in 1:4
        vertex_a = simplex[a + 1]
        vertex_b = simplex[b + 1]
        radial_a = edge_sq(edge_map, reference, vertex_a)
        radial_b = edge_sq(edge_map, reference, vertex_b)
        if a == b
            opposite = 0.0
        else
            opposite = edge_sq(edge_map, vertex_a, vertex_b)
        end
        # G_ab = (s_0a + s_0b - s_ab)/2.
        gram[a, b] = 0.5 * (radial_a + radial_b - opposite)
    end
    return gram
end

# metric reconstruction

function reconstruct_metric_from_gram(gram::AbstractMatrix{<:Real}, edge_matrix::AbstractMatrix{<:Real})
    if size(gram, 1) != size(gram, 2)
        error("Square matrix required")
    end
    
    if size(edge_matrix, 2) != size(gram, 1)
        error("Same dimensions required")
    end
    
    inv_edge_matrix = inv(Matrix{Float64}(edge_matrix))
    return transpose(inv_edge_matrix) * Matrix{Float64}(gram) * inv_edge_matrix
end

