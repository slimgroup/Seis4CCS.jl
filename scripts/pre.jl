using Pkg
Pkg.rm("AzureClusterlessHPC")
Pkg.develop(url="https://github.com/microsoft/AzureClusterlessHPC.jl.git")
Pkg.rm("JUDI4Cloud")
Pkg.develop(url="https://github.com/slimgroup/JUDI4Cloud.jl.git")