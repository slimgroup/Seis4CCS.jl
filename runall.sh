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
nsrc=${nsrc:=8}
vm=${vm:=Standard_F4}
nth=${nth:=4}
niter=${niter:=16}

julia ConcToV.jl --nv $L
julia GenData.jl --nv $L --nsrc $nsrc --vm $vm --nth $nth
julia JRM.jl --nv $L --vm $vm --nth $nth --niter $niter
