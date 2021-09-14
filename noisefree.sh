cd scripts

L=64
nsrc=16
niter=8
bs=4
vm=Standard_E8s_v3
nth=2

julia JRMsim.jl --nv $L --nsrc $nsrc --vm $vm --nth $nth --niter $niter --bs $bs