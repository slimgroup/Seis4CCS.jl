function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--nv"
            help = "Number of vintages"
            arg_type = Int
            default = 2
        "--nsrc"
            help = "Number of sources per vintage"
            arg_type = Int
            default = 8
        "--vm"
            help = "Type of virtual machine"
            arg_type = String
            default = "Standard_F4"
        "--nth"
            help = "Number of threads on a single node"
            arg_type = Int
            default = 4
        "--niter"
            help = "JRM iterations"
            arg_type = Int
            default = 8
        "--bs"
            help = "batchsize in JRM iterations"
            arg_type = Int
            default = 4
        "--snr"
            help = "SNR of noisy data"
            arg_type = Float64
            default = 0.0
    end
    return parse_args(s)
end