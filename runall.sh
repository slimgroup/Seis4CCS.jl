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

julia ConcToV.jl --nv $L
julia GenData.jl --nv $L
julia JRM.jl --nv $L
