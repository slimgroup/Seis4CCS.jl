using JUDI, FFTW, JLD2, Random, LinearAlgebra, Statistics

Random.seed!(1234)

using ArgParse
include("../utils/parse_cmd.jl")
parsed_args = parse_commandline()
L = parsed_args["nv"]
nsrc = parsed_args["nsrc"]
snr = parsed_args["snr"]

JLD2.@load "../data/dobs$(L)vint$(nsrc)nsrc.jld2" dobs_stack q_stack

scale = 10^(-snr/20)

noise_stack = deepcopy(dobs_stack)

for i = 1:length(noise_stack)
    for j = 1:noise_stack[i].nsrc
        e = randn(Float32,size(dobs_stack[i].data[j]))
        noise_stack[i].data[j] = real.(ifft(fft(e).*fft(q_stack[i].data[j])))
        noise_stack[i] = noise_stack[i]/norm(noise_stack[i])*scale*norm(dobs_stack[i])
    end
end

dobs_stack = dobs_stack + noise_stack

JLD2.@load "../data/dobs$(L)vint$(nsrc)nsrc.jld2" dobs_stack q_stack