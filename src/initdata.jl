# initdata.jl
# Initialization of the initial hyersurface slice
# two methods: one for a test gauge wave, which is fairly simple with coordinate rewrite
# and a general one for any metric which will prbably use the geodesic equations and path or something equivilent
#
#  

##########################################
# Gauge wave code
##########################################

const A = 0.1
const Nx = 20
const Ny = 20
const Nz = 20

const Lx = 1.0
const Ly = 1.0
const Lz = 1.0

const nsub=200

wrap_index(i::Int, n::Int)=mod(i,n)

function vertex_id(i::Int, j::Int, k::Int) # might make more versions of this, kind of doesn't really matter how you do it as long as it is unique
    ii=wrap_index(i,Nx);
    jj=wrap_index(j,Ny);
    kk=wrap_index(k,Nz);


    return 1 + ii + Nx*(jj + Ny*kk)
    
end

function build_cartesian_vertex_coords()
    
    dx=Lx/Nx;
    dy=Ly/Ny;
    dz=Lz/Nz;
    
     coords = Vector{SVector{3,Float64}}(undef, Nx * Ny * Nz)

    for k in 0:(Nz - 1), j in 0:(Ny - 1), i in 0:(Nx - 1)
        vid = vertex_id(i, j, k)
        coords[vid] = @SVector [i * dx, j * dy, k * dz]
    end

    return coords
    
end

function gauge_wave_H(x::Float64)
    return 1.0 + A * sin(2.0 * pi * x / Lx)
end

# we define the coordinate X = integral_0^x sqrt(H(ξ)) dξ so that dX^2+dy^2+dz^2 is the metric; in this way we can easily do the initialization
# here we integrate X using Simpson's rule
# also julia is pretty cool with the latex symbols 
function gauge_wave_X(x::Float64)
    if x == 0.0
        return 0.0
    elseif x < 0.0
        return -gauge_wave_X(-x)
    end

    h = x / nsub
    s = sqrt(gauge_wave_H(0.0)) + sqrt(gauge_wave_H(x))

    for m in 1:(nsub - 1)
        ξ = m * h
        coeff = isodd(m) ? 4.0 : 2.0
        s += coeff * sqrt(gauge_wave_H(ξ))
    end

    return (h / 3.0) * s
end

# since we will be building the same edge lengths at the same X location, we should store them somewhere
function build_gauge_wave_X_table()
    dx = Lx / Nx
    Xvals = Vector{Float64}(undef, Nx)

    for i in 0:(Nx - 1)
        x = i * dx
        Xvals[i + 1] = gauge_wave_X(x)
    end

    CX = gauge_wave_X(Lx) # total X distance, used for wrappin
    return Xvals, CX
end


function wrapped_delta_1d(x1::Float64, x2::Float64, L::Float64)

    d = x2-x1
    
    if d< -L*0.5
        d += L
    elseif d >= L*0.5
        d -= L
    end
    
    return d
end

# displacement in intrinsic X, using precomputed X values
function wrapped_delta_X(
    x1::Float64, x2::Float64,
    Xvals::Vector{Float64}, CX::Float64
)
    
    dx = Lx / Nx

    i1 = round(Int, x1 / dx)
    i2 = round(Int, x2 / dx)

    i1 = wrap_index(i1, Nx)
    i2 = wrap_index(i2, Nx)

    X1 = Xvals[i1 + 1]
    X2 = Xvals[i2 + 1]
    
    return wrapped_delta_1d(X1,X2,CX)
end

# wrapped total displacement with the X coordinate
function wrapped_intrinsic_displacement(
    p::SVector{3,Float64},
    q::SVector{3,Float64},
    Xvals::Vector{Float64}, CX::Float64
)
    dX = wrapped_delta_X(p[1], q[1], Xvals, CX)
    dy = wrapped_delta_1d(p[2], q[2], Ly)
    dz = wrapped_delta_1d(p[3], q[3], Lz)
    return @SVector [dX, dy, dz]
end






# intrinsic squared edge length for the gauge-wave slice using X
function edge_sq_from_gauge_wave_integral(
    p::SVector{3,Float64},
    q::SVector{3,Float64},
    Xvals::Vector{Float64}, CX::Float64
)
    d = wrapped_intrinsic_displacement(p, q, Xvals, CX)
    return dot(d, d)
end


# tetrahedral grid generation
# return the 8 vertices of one cell
function cell_vertices(i::Int, j::Int, k::Int)
    v000 = vertex_id(i,     j,     k)
    v100 = vertex_id(i + 1, j,     k)
    v010 = vertex_id(i,     j + 1, k)
    v110 = vertex_id(i + 1, j + 1, k)

    v001 = vertex_id(i,     j,     k + 1)
    v101 = vertex_id(i + 1, j,     k + 1)
    v011 = vertex_id(i,     j + 1, k + 1)
    v111 = vertex_id(i + 1, j + 1, k + 1)

    return (v000, v100, v010, v110, v001, v101, v011, v111)
end

# 6-tetra decomposition of a cube
function cube_to_tetras(
    v000::Int, v100::Int, v010::Int, v110::Int,
    v001::Int, v101::Int, v011::Int, v111::Int
)
    return [
        tetrakey(v000, v100, v110, v111),
        tetrakey(v000, v100, v101, v111),
        tetrakey(v000, v001, v101, v111),
        tetrakey(v000, v001, v011, v111),
        tetrakey(v000, v010, v011, v111),
        tetrakey(v000, v010, v110, v111),
    ]
end

# all tetrahedra and all unique edges
function build_periodic_tetrahedral_grid()
    tetra_set = Set{NTuple{4,Int}}()

    for k in 0:(Nz - 1), j in 0:(Ny - 1), i in 0:(Nx - 1)
        verts = cell_vertices(i, j, k)
        for tet in cube_to_tetras(verts...)
            push!(tetra_set, tet)
        end
    end

    tetras = collect(tetra_set)
    sort!(tetras)

    edge_set = Set{NTuple{2,Int}}()
    for tet in tetras
        a, b, c, d = tet
        push!(edge_set, edgekey(a, b))
        push!(edge_set, edgekey(a, c))
        push!(edge_set, edgekey(a, d))
        push!(edge_set, edgekey(b, c))
        push!(edge_set, edgekey(b, d))
        push!(edge_set, edgekey(c, d))
    end

    edges = collect(edge_set)
    sort!(edges)

    num_vert = Nx * Ny * Nz
    return num_vert, edges, tetras
end

# Build the fixed slice template and coordinates for the gauge-wave test.
function build_gauge_wave_template()

    coords = build_cartesian_vertex_coords()
    num_vert, edges, tetras = build_periodic_tetrahedral_grid()
    template = build_template(num_vert, edges, tetras)
    return template, coords
end

# Initialize the current slice squared edge lengths from the gauge-wave metric.
function initialize_gauge_wave_edge_s(
    template::SliceTemplate,
    vertex_coords::Vector{SVector{3,Float64}}
)

    Xvals, CX = build_gauge_wave_X_table()

    edge_s = Vector{Float64}(undef, length(template.edges))

    for (eid, (a, b)) in enumerate(template.edges)
        p = vertex_coords[a]
        q = vertex_coords[b]
        edge_s[eid] = edge_sq_from_gauge_wave_integral(p, q, Xvals, CX)
    end

    return edge_s
end

# One-call helper: build template, coordinates, current edge lengths, and state.
function build_gauge_wave_initial_state(
    slice_index::Int = 0
 )
   
    template, coords = build_gauge_wave_template()
    current_edge_s = initialize_gauge_wave_edge_s(template, coords)
    state = initialize_state(template, current_edge_s; slice_index = slice_index)

    return template, coords, state
end

















