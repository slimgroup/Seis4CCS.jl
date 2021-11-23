
using Random, Images, JLD2, LinearAlgebra
using JOLI, Statistics, FFTW
using JUDI
using Polynomials
using ArgParse
include("../utils/parse_cmd.jl")
parsed_args = parse_commandline()
L = parsed_args["nv"]

JLD2.@load "../models/Compass_tti_625m.jld2"
JLD2.@load "../data/Conc.jld2" conc

extentx = (n[1]-1)*d[1]
extentz = (n[2]-1)*d[2]

v = 1f3.*sqrt.(1f0./m)

idx_wb = find_water_bottom(rho.-rho[1,1])
idx_ucfmt = find_water_bottom((v.-3500f0).*(v.>3500f0))

Kh = zeros(Float32, n)
for i = 1:n[1]
    Kh[i,1:idx_wb[i]-1] .= 1f-4  # water layer
    Kh[i,idx_wb[i]:idx_ucfmt[i]-1] .= 1.03*1f-3*v[i,idx_wb[i]:idx_ucfmt[i]-1].+15f0
    Kh[i,idx_ucfmt[i]:idx_ucfmt[i]+7] .= 1f-4
    Kh[i,idx_ucfmt[i]+8:end] .= 1.03*1f-3*v[i,idx_ucfmt[i]+8:end].+200f0
end

phi = zeros(Float32,n)
for i = 1:n[1]
    for j = 1:n[2]
        p = Polynomial([-0.0314^2*Kh[i,j],2*0.0314^2*Kh[i,j],-0.0314^2*Kh[i,j],1.527^2])
        phi[i,j] = minimum(real(roots(p)[findall(real(roots(p)).== roots(p))]))
    end
    for j = idx_ucfmt[i]:idx_ucfmt[i]+7
        phi[i,idx_ucfmt[i]:idx_ucfmt[i]+7] = Float32.(range(0.056,stop=0.1,length=8))
    end
end

vp = convert(Array{Float32},v)
vs = vp ./ sqrt(3f0)

idx_wb = maximum(find_water_bottom(rho.-rho[1,1]))

function Patchy(sw, vp, vs, rho, phi; bulk_min = 36.6f9, bulk_fl1 = 2.735f9, bulk_fl2 = 0.125f9, ρw = 7.766f2, ρo = 1.053f3)
    bulk_sat1 = rho .* (vp.^2f0 - 4f0/3f0 .* vs.^2f0) * 1f3
    shear_sat1 = rho .* (vs.^2f0) * 1f3

    bulk_min = zeros(Float32,size(bulk_sat1))

    bulk_min[findall(v.>=3350)] .= 50f9   # mineral bulk moduli
    bulk_min[findall(v.<3350)] .= 1.2 * bulk_sat1[findall(v.<3350)] # mineral bulk moduli

    patch_temp = bulk_sat1 ./(bulk_min .- bulk_sat1) - bulk_fl1 ./ phi ./ (bulk_min .- bulk_fl1) + bulk_fl2 ./ phi ./ (bulk_min .- bulk_fl2)

    for i = 1:n[1]
        for j = idx_ucfmt[i]:idx_ucfmt[i]+8
            patch_temp[i,idx_ucfmt[i]:idx_ucfmt[i]+8] = patch_temp[i,idx_ucfmt[i]+9:idx_ucfmt[i]+17]
        end
    end
    
    bulk_sat2 = bulk_min./(1f0./patch_temp .+ 1f0)
    bulk_sat2[findall(bulk_sat2-bulk_sat1.>0)] = bulk_sat1[findall(bulk_sat2-bulk_sat1.>0)]

    bulk_new = (bulk_sat1+4f0/3f0*shear_sat1).*(bulk_sat2+4f0/3f0*shear_sat1) ./( (1f0.-sw).*(bulk_sat2+4f0/3f0*shear_sat1) 
    + sw.*(bulk_sat1+4f0/3f0*shear_sat1) ) - 4f0/3f0*shear_sat1

	bulk_new[:,1:idx_wb] = bulk_sat1[:,1:idx_wb]
    bulk_new[findall(sw.==0)] = bulk_sat1[findall(sw.==0)]
    rho_new = rho + phi .* sw * (ρw - ρo) / 1f3
    rho_new[findall(sw.==0)] = rho[findall(sw.==0)]
    Vp_new = sqrt.((bulk_new+4f0/3f0*shear_sat1)./rho_new/1f3)
    Vp_new[findall(sw.==0)] = vp[findall(sw.==0)]
    Vs_new = sqrt.((shear_sat1)./rho_new/1f3)

    return Vp_new, Vs_new, rho_new
end

indices = Int.(round.(range(1, stop=size(conc,1), length=L)))
sw_stack = [Float32.(imresize([zeros(Float32,19,325);conc[indices[i],:,:]]',n)) for i = 1:L]

vp_stack = [Patchy(sw_stack[i], vp, vs, rho, phi)[1] for i = 1:L]
rho_stack = [Patchy(sw_stack[i], vp, vs, rho, phi)[3] for i = 1:L]

JLD2.@save "../models/timelapsevrho$(L)vint_finer.jld2" vp_stack rho_stack
