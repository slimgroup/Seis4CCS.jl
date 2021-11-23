
using PyPlot
using Random, Images, JLD2, LinearAlgebra
using JOLI, Statistics, FFTW
using Printf
using JUDI4Cloud
using MECurvelets

using ArgParse
include("../utils/parse_cmd.jl")
parsed_args = parse_commandline()
nv = 16
vm = "Standard_E16s_v3"
nth = 2
batchsize = 8
nsrc = 16

using JUDI4Cloud
creds=joinpath(pwd(),"..","credentials.json")
Random.seed!(1234);

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
JLD2.@load "../data/dobs$(nv)vint$(nsrc)nsrcnoisefree.jld2" dobs_stack q_stack

C = joMECurvelet2D(n; DDT=Float32, RDT=Complex{Float32})
support = BitVector(undef, size(C,1))   # initially zero

idx_wb = find_water_bottom(rho.-rho[1,1])

v = 1f0./sqrt.(m)
v0 = deepcopy(v)
v0[:,maximum(idx_wb)+1:end] = 1f0./convert(Array{Float32,2},imfilter(1f0./v[:,maximum(idx_wb)+1:end], Kernel.gaussian(10)))
m0 = 1f0./v0.^2f0

rho0 = deepcopy(rho)
rho0[:,maximum(idx_wb)+1:end] = 1f0./convert(Array{Float32,2},imfilter(1f0./rho[:,maximum(idx_wb)+1:end], Kernel.gaussian(10)))

model0 = Model(n,d,o,m0; rho=rho0, nb=80)

opt = JUDI.Options(isic=true)

# Preconditioner

Tm = judiTopmute(model0.n, maximum(idx_wb), 3)  # Mute water column
S = judiDepthScaling(model0)
Mr = Tm*S

niter = 8

for L = 2:nv

    init_culsterless(L; credentials=creds, vm_size=vm, pool_name="JRM$(L)vintages", verbose=1, nthreads=nth, auto_scale=false, n_julia_per_instance=batchsize)
    @printf("Simultaneous source JRM Iteration at vintage %d \n", L)
    flush(Base.stdout)

    w = ones(Float32, length(support)) - 0.5f0 .* support # support is 0.5, non-support is 1
    W = joDiag(1f0./w; DDT=Complex{Float32}, RDT=Complex{Float32})
    x = [zeros(Complex{Float32}, size(C,1)) for i = 1:L+1]
    z = [zeros(Complex{Float32}, size(C,1)) for i = 1:L+1]
    γ   = Float32.(sqrt(L)) # hyperparameter to tune

    # sim src acquisition
    xsrc_stack = [[[q_stack[i].geometry.xloc[s][1] for s = 1:dobs_stack[i].nsrc] for k = 1:batchsize] for i = 1:L]
    ysrc_stack = [[[0.0f0] for k = 1:batchsize] for i = 1:L]
    zsrc_stack = [[[q_stack[i].geometry.zloc[s][1] for s = 1:dobs_stack[i].nsrc] for k = 1:batchsize] for i = 1:L]

    src_geometry_stack = [Geometry(xsrc_stack[i], ysrc_stack[i], zsrc_stack[i]; dt=q_stack[i].geometry.dt[1], t=q_stack[i].geometry.t[1]) for i = 1:L]

    # Main loop

    for  j=1:niter

        # Set up weights for current iteration
        weights = [randn(Float32, batchsize, nsrc) for i = 1:L]

        # Create wavelet
        wavelet = [[hcat(q_stack[i].data...) .* weights[i][k:k,:] for k = 1:batchsize] for i = 1:L]

        # Create sim src
        q_j = [judiVector(src_geometry_stack[i], wavelet[i]) for i = 1:L]  # super shot for current iteration

        # Create sim data
        data_j = [[sum(weights[i][k,:].*dobs_stack[i].data) for k = 1:batchsize] for i = 1:L]
        dobs_j = [judiVector(Geometry(dobs_stack[i].geometry.xloc[1],dobs_stack[i].geometry.yloc[1],dobs_stack[i].geometry.zloc[1]; dt=dobs_stack[i].geometry.dt[1], t=dobs_stack[i].geometry.t[1], nsrc=batchsize),data_j[i]) for i = 1:L]

        dmx = [Mr*C'*W*(1f0/γ*x[1]+x[i+1]) for i = 1:L]
    
        phi, g1 = lsrtm_objective([model0 for i = 1:L], q_j, dobs_j, dmx; options=opt, nlind=true)
        g = [1f0/γ*W*C*Mr'*vec(sum(g1)), [W*C*Mr'*vec(g1[i]) for i=1:L]...]

	    @printf("At iteration %d function value is %2.2e \n", j, phi)
	    flush(Base.stdout)
	    
        # Step size and update variable
        t = 2f-3*phi/norm(g)^2f0; # fixed step
        z -= t .* g;
        if L == 2 && j == 1
	        global threshold = 30*[Float32.(quantile(abs.(z[i]), .98)) for i = 1:L+1]
            global threshold1 = threshold[1]
            global threshold2 = maximum(threshold[2:end])
        end
        if j == 1
            z -= 5 * t .* g;
        end
        x[1] = soft_thresholding(z[1], threshold1)
        x[2:L+1] = [soft_thresholding(z[i], threshold2) for i = 2:L+1]
    end
    finalize_culsterless()
    JLD2.@save "../results/JRMsimnv$(nv)no$(L)IncreJRMwl1.jld2" x z threshold support
    for i = 1:L+1
        global support = support .| (abs.(x[i]).>0f0)
    end
end
