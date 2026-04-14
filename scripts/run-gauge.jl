#push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using ReggePIES


template, coords, state = ReggePIES.build_gauge_wave_initial_state()

println("some tetras:")
for t in template.tetras[1:150]
    println(t)
end

println("\n some edges:")
for e in template.edges[1:150]
    println(e)
end

println("\n some edge lengths:")
for i in 1:150
    println(template.edges[i], "  s = ", state.current.edge_s[i])
end
