
using PyPlot
using Random, Images, JLD2, LinearAlgebra
using JOLI, Statistics, FFTW
using Printf
using JUDI4Cloud

using ArgParse
include("../utils/parse_cmd.jl")
parsed_args = parse_commandline()
L = parsed_args["nv"]
vm = parsed_args["vm"]
niter = parsed_args["niter"]
nth = parsed_args["nth"]
nsrc = parsed_args["nsrc"]
batchsize = parsed_args["bs"]

Random.seed!(1234);

creds=joinpath(pwd(),"..","credentials.json")
init_culsterless(batchsize*L; credentials=creds, vm_size=vm, pool_name="IndpRecovery", verbose=1, nthreads=nth, auto_scale=false)

JLD2.@load "../models/Compass_tti_625m.jld2"
JLD2.@load "../models/timelapsevrho$(L)vint.jld2" vp_stack rho_stack
JLD2.@load "../data/dobs$(L)vint$(nsrc)nsrc.jld2" dobs_stack q_stack

idx_wb = find_water_bottom(rho.-rho[1,1])

nsrc = dobs_stack[1].nsrc

v0_stack = deepcopy(vp_stack./1f3)
for i = 1:L
    v0_stack[i][:,maximum(idx_wb)+1:end] .= 1f0./convert(Array{Float32,2},imfilter(1f3./vp_stack[i][:,maximum(idx_wb)+1:end], Kernel.gaussian(10)))
end
m0_stack = [1f0./v0_stack[i].^2f0 for i = 1:L]

model0_stack = [Model(n,d,o,m0_stack[i]; nb=80) for i = 1:L]

opt = JUDI.Options(isic=true)

# Preconditioner

Tm = judiTopmute(model0_stack[1].n, maximum(idx_wb), 3)  # Mute water column
S = judiDepthScaling(model0_stack[1])
Mr = Tm*S

# Linearized Bregman parameters
nn = prod(model0_stack[1].n)
x = zeros(Float32, nn)
z = zeros(Float32, nn)

fval = zeros(Float32, niter)

# Soft thresholding functions and Curvelet transform
soft_thresholding(x::Array{Float64}, lambda) = sign.(x) .* max.(abs.(x) .- convert(Float64, lambda), 0.0)
soft_thresholding(x::Array{Float32}, lambda) = sign.(x) .* max.(abs.(x) .- convert(Float32, lambda), 0f0)

C0 = joCurvelet2D(2*n[1], 2*n[2]; zero_finest = false, DDT = Float32, RDT = Float64)

function C_fwd(im, C, n)
	im = hcat(reshape(im, n), reshape(im, n)[:, end:-1:1])
    im = vcat(im, im[end:-1:1,:])
	coeffs = C*vec(im)/2f0
	return coeffs
end

function C_adj(coeffs, C, n)
	im = reshape(C'*coeffs, 2*n[1], 2*n[2])
	return vec(im[1:n[1], 1:n[2]] + im[1:n[1], end:-1:n[2]+1] + im[end:-1:n[1]+1, 1:n[2]] + im[end:-1:n[1]+1, end:-1:n[2]+1])/2f0
end

C = joLinearFunctionFwd_T(size(C0, 1), n[1]*n[2],
                          x -> C_fwd(x, C0, n),
                          b -> C_adj(b, C0, n),
                          Float32,Float64, name="Cmirrorext")

src_list = [collect(1:nsrc) for i = 1:L]

ps = 0
γ  = 1f0 # hyperparameter to tune

x = [zeros(Float32, size(C,2)) for i=1:L];
z = [zeros(Float64, size(C,1)) for i=1:L];

tau = [zeros(Float32,size(C,1)) for i=1:L];
flag = [BitArray(undef,size(C,1)) for i=1:L];
sumsign = [zeros(Float32,size(C,1)) for i=1:L]

lambda = zeros(Float64,L)

# Main loop

for  j=1:niter

	# Main loop			   
    @printf("Independent Recovery Iteration: %d \n", j)
    flush(Base.stdout)
				
    global inds = [zeros(Int, batchsize) for i=1:L]

    for i = 1:L
        length(src_list[i]) < batchsize && (global src_list[i] = collect(1:nsrc))
		src_list[i] = src_list[i][randperm(MersenneTwister(i+2000*j),length(src_list[i]))]
		global inds[i] = [pop!(src_list[i]) for b=1:batchsize]
		println("Vintage $(i) Imaging source $(inds[i])")
    end

    source = [q_stack[i][inds[i]] for i = 1:L]
    dObs = [dobs_stack[i][inds[i]] for i = 1:L]
    dmx = [Mr*x[i] for i = 1:L]
    
    phi, g1 = lsrtm_objective(model0_stack, source, dObs, dmx; options=opt, nlind=false)
    g = [C*Mr'*vec(g1[i]) for i=1:L]

	@printf("At iteration %d function value is %2.2e \n", j, phi)
	flush(Base.stdout)
	# Step size and update variable

	global t = 2f-5/12.5 # fixed step

    # anti-chatter
	for i = 1:L
		global sumsign[i] = sumsign[i] + sign.(g[i])
		global tau[i] .= t
		global tau[i][findall(flag[i])] = deepcopy((t*abs.(sumsign[i])/j)[findall(flag[i])])
		global z[i] -= tau[i] .* g[i]
	end

	(j==1) && global lambda = [quantile(abs.(vec(z[1])), .8) for i = 1:L]	# estimate thresholding parameter at 1st iteration

	# Update variables and save snapshot
	global x = [adjoint(C)*soft_thresholding(z[i], lambda[i]) for i = 1:L]
	for i = 1:L
		global flag[i] = flag[i] .| (abs.(z[i]).>=lambda[i])     # check if ever pass the threshold
	end
    JLD2.@save "../results/IndpRecov$(j)Iter$(L)vintages$(nsrc)nsrc.jld2" x z g lambda phi
end
