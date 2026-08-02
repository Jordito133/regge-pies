using ReggePIES
using Printf

amp = 0.1
wave = 1.0
number_frames = 300

model = ReggePIES.gauge_wave_metric(
	amplitude=amp,
	wavelength=wave,
)

grid = ReggePIES.GridSpec(
	(300, 300, 300),
	(1.2, 1.2, 1.2),
	(true, true, true),
)

frame_directory = joinpath(@__DIR__, "..", "output", "gauge-frames")
mkpath(frame_directory)

for frame in 0:(number_frames - 1)
    time = wave * frame / number_frames

    lapse_values = ReggePIES.lapse_plane(model, grid, time)
    image = ReggePIES.grayscale_grid(
        lapse_values;
        cell_pixels=2,
    )

    filename = @sprintf("frame_%05d.pgm", frame + 1)
    path = joinpath(frame_directory, filename)
    ReggePIES.write_pgm(path, image)

    println("Wrote $filename")
end