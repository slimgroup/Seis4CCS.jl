#!/bin/bash

cd scripts

L=${L:=64}
nsrc=${nsrc:=16}
vm=${vm:=Standard_E8s_v3}
nth=${nth:=8}
niter=${niter:=6}
bs=${bs:=8}
snr=${snr:=0}
gamma=${gamma:=16}

julia -p 2 -L ~/startup.jl GenNonlinearData.jl --nv $L --nsrc $nsrc --vm $vm --nth $nth
julia -p 2 -L ~/startup.jl JRMsim1bg.jl --nv $L --nsrc $nsrc --vm $vm --nth $nth --niter $niter --bs $bs --gamma $gamma