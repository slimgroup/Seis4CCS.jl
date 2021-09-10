
function plot_image(image, spacing; perc=98, cmap="Greys", o=(0, 0),
    interp="hanning", aspect=nothing, d_scale=1.5,
    name=nothing, units="km", new_fig=true, save=nothing, a=nothing,start_grid=nothing,alpha=1,labelsize=12)
nx, nz = size(image)
dx, dz = spacing
ox, oz = o
depth = range(oz, oz + (nz - 1)*spacing[2], length=nz).^d_scale

scaled = image .* depth'
if isnothing(a)
    a = quantile(abs.(vec(scaled)), perc/100)
end
start_x = start_z = 1
if !isnothing(start_grid)
    start_x = start_grid[1]
    start_z = start_grid[2]
end
extent = [(start_x-1)*dx/1f3+ox/1f3, (start_x-1)*dx/1f3+(ox+ (nx-1)*dx)/1f3, (start_z-1)*dz/1f3+(oz+(nz-1)*dz)/1f3, (start_z-1)*dz/1f3+oz/1f3]
isnothing(aspect) && (aspect = .5 * nx/nz)

if new_fig
figure()
end
PyPlot.rc("font", family="serif"); PyPlot.rc("xtick", labelsize=labelsize); PyPlot.rc("ytick", labelsize=labelsize)
imshow(scaled', vmin=-a, vmax=a, cmap=cmap, aspect=aspect,
interpolation=interp, extent=extent,alpha=alpha)
xlabel("X [$units]", fontsize=labelsize)
ylabel("Z [$units]", fontsize=labelsize)
if !isnothing(name)
    title("$name")
end

if ~isnothing(save)
save == True ? filename=name : filename=save
savefig(filename, bbox_inches="tight", dpi=150)
end
return a
end


function plot_data(data, dt; perc=98, cmap="seismic",
    interp="hanning", aspect=nothing,
    name="shot record", units="s", new_fig=true, save=nothing, a=nothing,labelsize=12)
    nt, nrec = size(data)

if isnothing(a)
    a = quantile(abs.(vec(data)), perc/100)
end
extent = [1, nrec, (nt-1)*dt/1f3,0]
isnothing(aspect) && (aspect = 2*nrec/((nt-1)*dt/1f3))

if new_fig
figure()
end
PyPlot.rc("font", family="serif"); PyPlot.rc("xtick", labelsize=labelsize); PyPlot.rc("ytick", labelsize=labelsize)
imshow(data, vmin=-a, vmax=a, cmap=cmap, aspect=aspect,
interpolation=interp, extent=extent)
xlabel("Receiver no.", fontsize=labelsize)
ylabel("Time [$(units)]", fontsize=labelsize)
title("$name")

if ~isnothing(save)
save == True ? filename=name : filename=save
savefig(filename, bbox_inches="tight", dpi=150)
end
return a
end
