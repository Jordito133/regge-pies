#push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using ReggePIES


num_vert = 6
edges = [
    (1,2),
    (1,3),
    (1,4),
    (2,3),
    (2,4),
    (2,5),
    (2,6),
    (3,4),
    (3,5),
    (3,6),
    (4,5),
    (4,6),
    (5,6),
    

  ]
tetras = [
    (1,2,3,4),
    (2,3,4,5),
    (2,4,5,6),
    (2,3,4,6),
    (3,4,5,6),
    
        
]

template = ReggePIES.build_template(num_vert, edges, tetras)

println("Neighbors:")
for v in 1:num_vert
    println("v = $v : ", template.vertex_neighbors[v])
end

println("\nVertex to edges:")
for v in 1:num_vert
    println("v = $v : ", template.vertex_to_edges[v])
end

println("\nVertex to tetras:")
for v in 1:num_vert
    println("v = $v : ", template.vertex_to_tetras[v])
end

println("\nEdge to triangles:")
for (eid, tri_list) in enumerate(template.edge_to_triangles)
    println("edge $eid $(template.edges[eid]) : $tri_list")
end

groups = ReggePIES.greedy_groups(template)
println("\nParallel groups:")
for (i,g) in enumerate(groups)
    println("group $i : $g")
end

 current_edge_s = [
    2.0,  
    3.0,  
    2.0,  
    2.0,  
    2.0,  
    2.0,  
    2.0,
    2.0,
    2.0,
    2.0,
    3.0,
    4.0,
    -1.0,
    
]

state = ReggePIES.initialize_state(template, current_edge_s; slice_index=0)

println("\nSlice edge lengths:")
for (eid, e) in enumerate(template.edges)
    println("edge $eid $e  s = ", state.current.edge_s[eid])
end
