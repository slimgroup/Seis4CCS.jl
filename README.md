# Seis4CCS.jl

This repository contains scripts to run seismic monitoring of CCS with joint recovery on Azure batch. These scripts are suggested to run from docker image, which will be provided upon request.

## TO-DO

First, please put your Azure credential file in this directory (`Seis4CCS.jl/credentials.json`). e.g. do

```bash
touch credentials.json
vi credentials.json
```

and copy-paste your credentials to this file.

Then, under `Seis4CCS.jl`, run

```bash
L=x nsrc=y bash runall.sh
```

where integer `x` is the number of vintages in the experiment (maximum is 879), and integer `y` is number of jittered sample seismic sources in each vintage (different for each vintage, suggested to be more than `8`). `x` and `y` willl be provided by user. The generated seismic data will be saved in `./data` folder and recovery results will be saved in `./results` folder.

There is `rclone` installed in the docker image. Feel free to use it to upload the results to your preferred online platforms.

## Authors

If any question, please contact

Ziyi (Francis) Yin, ziyi.yin@gatech.edu    
Mathias Louboutin, mlouboutin3@gatech.edu    
Felix J. Herrmann, felix.herrmann@gatech.edu     