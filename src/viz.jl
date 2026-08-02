# viz.jl
# vizualization of the metrix in a 2D grid using the lapse
#
#

function lapse_plane(model::SpacetimeMetricModel, grid::GridSpec, time::Real; z_index::Int = cld(grid.N[3],2)) #z_index = grid.Nz/2
    
    values = Matrix{Float64}(undef, grid.N[2], grid.N[1])
    z=(z_index-1)*grid_spacing(grid, 3)
    for j in 0:(grid.N[2]-1)
        for i in 0: (grid.N[1]-1)
            point = SVector{3, Float64}(
                i*grid_spacing(grid,1),
                j*grid_spacing(grid, 2),
                z,
            )
            values[grid.N[2]-j , i+1] = lapse(model, time, point)
            #values on the screen are in an inverted coordinate frame
        end
    end
    return values
end


function grayscale_grid(values::AbstractMatrix{Float64};
    cell_pixels::Int = 8
)
    # will change these in code
    #########################
    min = 0.94
    max = 1.06
    ####################
    
    
    
    
    rows, columns = size(values)
    height = rows*cell_pixels
    width = columns*cell_pixels
    
    image = zeros(UInt8, height, width)
    
    for r in 1:rows
        for c in 1:columns
            norm=clamp((values[r, c]-min)/(max-min), 0.0, 1.0)
            shade = UInt8(round(Int, 255.0*norm ))
            frow = (r-1)*cell_pixels+1
            fcol = (c-1)*cell_pixels+1
            image[
            frow:(frow+cell_pixels-1),
            fcol:(fcol+cell_pixels-1),
            ] .= shade # all pixels in this area are the same color
        end
    end
return image
    
    
    
end

# we need to write the array of colors to an image - seems like PGM works

function write_pgm(path::AbstractString, image::AbstractMatrix{UInt8})
    mkpath(dirname(path))
    height, width = size(image)
    open(path, "w") do stream
        write(stream, "P5\n$width $height\n255\n")
   
        write(stream, permutedims(Matrix(image)))
    end
    return String(path)
end
