
using PyPlot
using Random, Images, JLD2, LinearAlgebra
using JOLI, Statistics, FFTW
using Printf
using JUDI
using MECurvelets

using ArgParse
include("../utils/parse_cmd.jl")
parsed_args = parse_commandline()
L = parsed_args["nv"]
vm = parsed_args["vm"]
niter = parsed_args["niter"]
nth = parsed_args["nth"]
nsrc = parsed_args["nsrc"]
batchsize = parsed_args["bs"]

import MECurvelets.soft_thresholding

function soft_thresholding(c::Vector{Complex{vDt}}, λ::vDt) where {vDt}
    out = exp.(angle.(c)im) .* max.(abs.(c) .- convert(vDt, λ), 0)
    return out
end

function soft_thresholding(c::Vector{Complex{vDt}}, λ::Vector{vDt}) where {vDt}
    out = exp.(angle.(c)im) .* max.(abs.(c) .- convert.(vDt, λ), 0)
    return out
end

Random.seed!(1234);

JLD2.@load "../models/Compass_tti_625m.jld2"
JLD2.@load "../models/timelapsevrho$(L)vint.jld2" vp_stack rho_stack
JLD2.@load "../data/dobs$(L)vint$(nsrc)nsrcnoisefree.jld2" dobs_stack q_stack

C = joMECurvelet2D(n; DDT=Float32, RDT=Complex{Float32})
support = BitVector(undef, size(C,1))   # initially zero

idx_wb = find_water_bottom(rho.-rho[1,1])

nsrc = dobs_stack[1].nsrc

v0_stack = deepcopy(vp_stack./1f3)
for i = 1:L
    v0_stack[i][:,maximum(idx_wb)+1:end] .= 1f0./convert(Array{Float32,2},imfilter(1f3./vp_stack[1][:,maximum(idx_wb)+1:end], Kernel.gaussian(10)))
end
m0_stack = [1f0./v0_stack[i].^2f0 for i = 1:L]

rho0_stack = deepcopy(rho_stack)
for i = 1:L
    rho0_stack[i][:,maximum(idx_wb)+1:end] .= 1f0./convert(Array{Float32,2},imfilter(1f0./rho_stack[1][:,maximum(idx_wb)+1:end], Kernel.gaussian(10)))
end

model0_stack = [Model(n,d,o,m0_stack[i];rho=rho0_stack[i], nb=80) for i = 1:L]

opt = JUDI.Options(isic=true)

# Preconditioner

Tm = judiTopmute(model0_stack[1].n, maximum(idx_wb), 3)  # Mute water column
S = judiDepthScaling(model0_stack[1])
Mr = Tm*S

# Soft thresholding functions and Curvelet transform
soft_thresholding(x::Array{Float64}, lambda) = sign.(x) .* max.(abs.(x) .- convert(Float64, lambda), 0.0)
soft_thresholding(x::Array{Float32}, lambda) = sign.(x) .* max.(abs.(x) .- convert(Float32, lambda), 0f0)

x = [zeros(Float32, size(C,2)) for i=1:L+1];
z = [zeros(Float64, size(C,1)) for i=1:L+1];

flag = [BitArray(undef,size(C,1)) for i=1:L+1];
sumsign = [zeros(Float32,size(C,1)) for i=1:L+1];

lambda = zeros(Float64,L+1);

# sim src acquisition
xsrc_stack = [[[q_stack[i].geometry.xloc[s][1] for s = 1:dobs_stack[i].nsrc] for k = 1:batchsize] for i = 1:L]
ysrc_stack = [[[0.0f0] for k = 1:batchsize] for i = 1:L]
zsrc_stack = [[[q_stack[i].geometry.zloc[s][1] for s = 1:dobs_stack[i].nsrc] for k = 1:batchsize] for i = 1:L]

src_geometry_stack = [Geometry(xsrc_stack[i], ysrc_stack[i], zsrc_stack[i]; dt=q_stack[i].geometry.dt[1], t=q_stack[i].geometry.t[1]) for i = 1:L]


# Main loop

for  j=1:niter

	# Main loop			   
    @printf("Simultaneous source JRM Iteration: %d \n", j)
    flush(Base.stdout)

    # Set up weights for current iteration
    weights = [randn(Float32, batchsize, nsrc) for i = 1:L]

    # Create wavelet
    wavelet = [[hcat(q_stack[i].data...) .* weights[i][k:k,:] for k = 1:batchsize] for i = 1:L]

    # Create sim src
    q_j = [judiVector(src_geometry_stack[i], wavelet[i]) for i = 1:L]  # super shot for current iteration

    # Create sim data
    data_j = [[sum(weights[i][k,:].*dobs_stack[i].data) for k = 1:batchsize] for i = 1:L]
    dobs_j = [judiVector(Geometry(dobs_stack[i].geometry.xloc[1],dobs_stack[i].geometry.yloc[1],dobs_stack[i].geometry.zloc[1]; dt=dobs_stack[i].geometry.dt[1], t=dobs_stack[i].geometry.t[1], nsrc=batchsize),data_j[i]) for i = 1:L]

    dmx = [Mr*(1f0/γ*x[1]+x[i+1]) for i = 1:L]
    
    phi1 = Array{Float32, 1}(undef, L)
    g1 = Array{Array{Float32,2}, 1}(undef, L)
    for i = 1:L
        phi1[i], g1[i] = lsrtm_objective(model0_stack[i], q_j[i], dobs_j[i], dmx[i]; options=opt, nlind=true)
    end
    phi = sum(phi1)
    g = [1f0/γ*C*Mr'*vec(sum(g1)), [C*Mr'*vec(g1[i]) for i=1:L]...]

	@printf("At iteration %d function value is %2.2e \n", j, phi)
	flush(Base.stdout)
	# Step size and update variable

	global t = γ^2f0*4f-7/(L+γ^2f0)*32/nsrc # fixed step

    # anti-chatter
    for i = 1:L+1
        global sumsign[i] = sumsign[i] + sign.(g[i])
        tau = t*ones(Float32,size(C,1))
        tau[findall(flag[i])] = deepcopy((t*abs.(sumsign[i])/j)[findall(flag[i])])
        global z[i] -= tau .* g[i]
    end

	(j==1) && global lambda = [quantile(abs.(vec(z[1])), .6), [quantile(abs.(vec(z[i])), .95) for i = 2:L+1]...]	# estimate thresholding parameter at 1st iteration
    lambda1 = maximum(lambda[2:end])
    for i = 2:L+1
        global lambda[i] = lambda1
    end

	# Update variables and save snapshot
	for i = 1:L+1
		global x[i] = adjoint(C)*soft_thresholding(z[i], lambda[i])
	end
    JLD2.@save "../results/JRMsim$(j)Iter$(L)vintages$(nsrc)nsrc.jld2" x z g lambda phi
end
