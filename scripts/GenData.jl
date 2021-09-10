using Random, Images, JLD2, LinearAlgebra
using JOLI, Statistics, FFTW
using JUDI4Cloud
using Printf

using ArgParse
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--nv"
            help = "Number of vintages"
            arg_type = Int
            default = 2
    end
    return parse_args(s)
end
parsed_args = parse_commandline()
L = parsed_args["nv"]

Random.seed!(1234)

JLD2.@load "../models/Compass_tti_625m.jld2"
JLD2.@load "../models/timelapsevrho$(L)vint.jld2" vp_stack rho_stack
idx_wb = find_water_bottom(rho.-rho[1,1])
include("../utils/Jitter.jl")

m_stack = [(1f3./vp_stack[i]).^2f0 for i = 1:L]

model = Model(n,d,o,m; rho=rho,nb=80)

dtS = dtR = 4f0
timeS = timeR = 2500f0

nsrc = 4
nrec = Int.(floor((n[1]-1)*d[1]))   # 1m

creds=joinpath(pwd(),"..","credentials.json")
init_culsterless(nsrc; credentials=creds, vm_size="Standard_F4", pool_name="Fq", verbose=1, nthreads=4)

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

ntComp = get_computational_nt(srcGeometry_stack[1], recGeometry, model)
info = Info(prod(n), nsrc, ntComp)

opt = JUDI.Options(free_surface=true)

F_stack = [judiProjection(info, recGeometry)*judiModeling(info, model; options=opt)*adjoint(judiProjection(info, srcGeometry_stack[i])) for i = 1:L]

dobs_stack = Array{judiVector{Float32,Array{Float32,2}},1}(undef, L)
for i = 1:L
    dobs_stack[i] = F_stack[i] * q_stack[i]
end

JLD2.@save "../data/dobs$(L)vint.jld2" dobs_stack q_stack

finalize_culsterless()
