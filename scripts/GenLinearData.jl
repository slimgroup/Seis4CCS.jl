using Random, Images, JLD2, LinearAlgebra
using JOLI, Statistics, FFTW
using JUDI4Cloud
using Printf
using Distributed

using ArgParse
include("../utils/parse_cmd.jl")
parsed_args = parse_commandline()
L = parsed_args["nv"]
nsrc = parsed_args["nsrc"]
vm = parsed_args["vm"]
nth = parsed_args["nth"]

creds=joinpath(pwd(),"..","credentials.json")
init_culsterless(nsrc*L; credentials=creds, vm_size=vm, pool_name="Jdm", verbose=1, nthreads=nth, auto_scale=false)

Random.seed!(1234)

JLD2.@load "../models/Compass_tti_625m.jld2"
JLD2.@load "../models/timelapsevrho$(L)vint.jld2" vp_stack rho_stack
idx_wb = find_water_bottom(rho.-rho[1,1])
include("../utils/Jitter.jl")

m_stack = [(1f3./vp_stack[i]).^2f0 for i = 1:L]

model_stack = [Model(n,d,o,m_stack[i]; nb=80) for i = 1:L]

dtS = dtR = 4f0
timeS = timeR = 2500f0

nrec = Int.(floor((n[1]-1)*d[1]))   # 1m

xsrc_stack = [convertToCell(ContJitter(n, d, nsrc)) for i = 1:L]
ysrc = convertToCell(range(0f0,stop=0f0,length=nsrc))
zsrc_stack = [convertToCell(range(10f0,stop=10f0,length=nsrc).+10f0*rand(Float32,nsrc).-5f0) for i = 1:L]

srcGeometry_stack = [Geometry(xsrc_stack[i], ysrc, zsrc_stack[i]; dt=dtS, t=timeS) for i = 1:L]
wavelet = ricker_wavelet(timeS,dtS,0.025f0)

q_stack = [judiVector(srcGeometry_stack[i], wavelet) for i = 1:L]

xrec = range(d[1],stop=(n[1]-1)*d[1],length=nrec)
yrec = 0f0
zrec = range((maximum(idx_wb)-1)*d[1]-2f0,stop=(maximum(idx_wb)-1)*d[1]-2f0,length=nrec)

recGeometry = Geometry(xrec,yrec,zrec; dt=dtR, t=timeR, nsrc=nsrc)

ntComp = get_computational_nt(srcGeometry_stack[1], recGeometry, model_stack[1])
info = Info(prod(n), nsrc, ntComp)

opt = JUDI.Options(isic=true)

F_stack = [judiProjection(info, recGeometry)*judiModeling(info, model_stack[i]; options=opt)*adjoint(judiProjection(info, srcGeometry_stack[i])) for i = 1:L]

v0_stack = deepcopy(vp_stack./1f3)
for i = 1:L
    v0_stack[i][:,maximum(idx_wb)+1:end] .= 1f0./convert(Array{Float32,2},imfilter(1f3./vp_stack[i][:,maximum(idx_wb)+1:end], Kernel.gaussian(10)))
end
m0_stack = [1f0./v0_stack[i].^2f0 for i = 1:L]

model0_stack = [Model(n,d,o,m0_stack[i]; nb=80) for i = 1:L]

F0_stack = [judiProjection(info, recGeometry)*judiModeling(info, model0_stack[i]; options=opt)*adjoint(judiProjection(info, srcGeometry_stack[i])) for i = 1:L]

J_stack = [judiJacobian(F0_stack[i], q_stack[i]) for i = 1:L]
dm_stack = [vec(m_stack[i]-m0_stack[i]) for i = 1:L]

dobs_stack = Array{judiVector{Float32,Array{Float32,2}},1}(undef, L)
@sync for i = 1:L
    @async dobs_stack[i] = J_stack[i] * dm_stack[i]
end

JLD2.@save "../data/dobs$(L)vint$(nsrc)nsrcnoisefree.jld2" dobs_stack q_stack
