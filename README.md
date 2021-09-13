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
L=a nsrc=b vm=c nth=d niter=e bs=f snr=g bash runall.sh
```

where integer `a` is the number of vintages in the experiment (maximum is 879), integer `b` is number of jittered sample seismic sources in each vintage (different for each vintage, suggested to be more than `8`), `c` is the type of virtual machine for seismic data simulation and imaging, `d` is the number of threads on each computing node in Azure batch, `e` is the number of joint recovery iterations, `f` is the batchsize in a JRM iteration, float `g` is the SNR of the noise to be added in the data. The default values of `a,b,c,d,e,f,g` are `2,8,Standard_F4,4,16,4,0.0`.

When the aformentioned program finishes, the generated seismic data will be saved in `./data` folder and each iteration of joint recovery results will be saved in `./results` folder. To test if everything runs without error in a small scale, you could do

```bash
bash runall.sh
```

Notice that if you set `L`, `nsrc`, `niter` to be very large, it may require massive disk space and memory on your local machine. It is also suggested to add `--oom-kill-disable` when you run the docker image.

There is `rclone` installed in the docker image. Feel free to use it to upload the results to your preferred online platforms.

## Authors

If any question, please contact

Ziyi (Francis) Yin, ziyi.yin@gatech.edu    
Mathias Louboutin, mlouboutin3@gatech.edu    
Felix J. Herrmann, felix.herrmann@gatech.edu     