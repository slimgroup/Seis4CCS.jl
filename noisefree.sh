#!/bin/bash

mkdir results

mkdir models
cd models
wget https://www.dropbox.com/s/n2cc5nntameh4rg/Compass_tti_625m.jld2
cd ..

mkdir data
cd data
wget https://www.dropbox.com/s/radan2hxgb5jnc5/Conc.jld2

cd ../scripts

L=${L:=2}
nsrc=${nsrc:=4}
vm=${vm:=Standard_E8s_v3}
nth=${nth:=8}
niter=${niter:=16}
bs=${bs:=4}
snr=${snr:=0}
gamma=${gamma:=1}

julia ConcToV.jl --nv $L
julia GenLinearData.jl --nv $L --nsrc $nsrc --vm $vm --nth $nth
julia JRMsim.jl --nv $L --nsrc $nsrc --vm $vm --nth $nth --niter $niter --bs $bs --gamma $gamma