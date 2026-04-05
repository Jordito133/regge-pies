# combinatorics.jl
# The combinatorics of the slices, used for all slices, fixed in place
# It should be the case that PIES does not change underlying topology
# Relations of edges, vertices, tetrahedra, triangles, and anything else that might be needed

#sorting edges, triangles and tetra from smallest to largest ids to remove repetitions 
edgekey(i::Int,j::Int) = i<j ? (i,j) : (j,i)
trianglekey(i::Int, j::Int, k::Int) = Tuple(sort((i, j, k)))
tetrakey(i::Int, j::Int, k::Int, l::Int) = Tuple(sort((i, j, k, l)))

# the main templete for a slice; internal topology and combinatorics stay unnefected after a PIES evolution step
struct SliceTemplate
    num_vert::Int
    edges::Vector{NTuple{2,Int}} #vertex ids
    tetras::Vector{NTuple{4,Int}} #vertex ids

    vertex_neighbors::Vector{Vector{Int}} #vector of vectors of vertex ids
    vertex_to_edges::Vector{Vector{Int}} # same for edge ids
    vertex_to_tetras::Vector{Vector{Int}} # tetras
    edge_to_triangles::Vector{Vector{NTuple{3,Int}}} # the triangles that an edge is a part of
end

# builder functions for init
function build_vertex_neighbors(num_vert::Int, edges::Vector{NTuple{2,Int}})
    neighbors = [Int[] for _ in 1:num_vert]
    for (a,b) in edges
        push!(neighbors[a], b)
        push!(neighbors[b], a)
    end
    return neighbors
end

function build_vertex_to_edges(num_vert::Int, edges::Vector{NTuple{2,Int}})
    v2e = [Int[] for _ in 1:num_vert]
    for (eid,(a,b)) in enumerate(edges)
        push!(v2e[a], eid)
        push!(v2e[b], eid)
    end
    return v2e
end

function build_vertex_to_tetras(num_vert::Int, tetras::Vector{NTuple{4,Int}})
    v2t = [Int[] for _ in 1:num_vert]
    for (tid, tet) in enumerate(tetras)
        for v in tet
            push!(v2t[v], tid)
        end
    end
    return v2t
end

function tetra_faces(t::NTuple{4,Int})
    a,b,c,d = t
    return [
        sort((a,b,c)),
        sort((a,b,d)),
        sort((a,c,d)),
        sort((b,c,d)),
    ]
end

function build_edge_to_triangles(edges::Vector{NTuple{2,Int}}, tetras::Vector{NTuple{4,Int}})
    edge_index = Dict{Tuple{Int,Int},Int}()
    for (eid,(a,b)) in enumerate(edges)
        edge_index[edgekey(a,b)] = eid
    end

    e2tri = [NTuple{3,Int}[] for _ in 1:length(edges)]

    seen = Set{NTuple{3,Int}}()
    for tet in tetras
        a,b,c,d = tet
        tris = [
            Tuple(sort((a,b,c))),
            Tuple(sort((a,b,d))),
            Tuple(sort((a,c,d))),
            Tuple(sort((b,c,d))),
        ]
        for tri in tris
            if tri in seen
                continue
            end
            push!(seen, tri)
            i,j,k = tri
            for e in [edgekey(i,j), edgekey(i,k), edgekey(j,k)]
                eid = edge_index[e]
                push!(e2tri[eid], tri)
            end
        end
    end

    return e2tri
end

function build_template(num_vert::Int,
                        edges::Vector{NTuple{2,Int}},
                        tetras::Vector{NTuple{4,Int}})
    neighbors = build_vertex_neighbors(num_vert, edges)
    v2e = build_vertex_to_edges(num_vert, edges)
    v2t = build_vertex_to_tetras(num_vert, tetras)
    e2tri = build_edge_to_triangles(edges, tetras)

    return SliceTemplate(num_vert, edges, tetras, neighbors, v2e, v2t, e2tri)
end
