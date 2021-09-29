#!/bin/bash

cd scripts

L=${L:=2}
nsrc=${nsrc:=64}
vm=${vm:=Standard_E8s_v3}
nth=${nth:=8}
niter=${niter:=8}
bs=${bs:=8}
snr=${snr:=0}
gamma=${gamma:=0.5}

julia -p 2 -L ~/startup.jl GenNonlinearData.jl --nv $L --nsrc $nsrc --vm $vm --nth $nth
julia -p 2 -L ~/startup.jl JRMsim1bg.jl --nv $L --nsrc $nsrc --vm $vm --nth $nth --niter $niter --bs $bs --gamma $gamma