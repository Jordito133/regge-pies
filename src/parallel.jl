# parallel.jl
# Building of the adjecency graphs and grouping vertecies in non-neighbouring groups 
# 

# Creating a graph for the slice
function slice_graph(template::SliceTemplate)
    g = SimpleGraph(template.num_vert)
    for (a,b) in template.edges
        add_edge!(g, a, b)
    end
    return g
end

# Greedy graph coloring/indexing
# NP complete problem, but this alg should be good enough i think
# 
function greedy_groups(template::SliceTemplate)
    g = slice_graph(template)
    color_of = Dict{Int,Int}() #color assigned to a given vertex

    for v in 1:template.num_vert
        used = Set{Int}()
        #collect the colors that are already used by neighbours of v
        for n in neighbors(g, v)
            if haskey(color_of, n)
                push!(used, color_of[n])
            end
        end
        # checking if all colors are used
        c = 1
        while c in used
            c += 1
        end
        color_of[v] = c
    end

    maxc = maximum(values(color_of))
    groups = [Int[] for _ in 1:maxc]
    for v in 1:template.num_vert
        push!(groups[color_of[v]], v)
    end

    return groups
end
