# Seis4CCS.jl

This repository contains scripts to run numerical simulations of seismic monitoring of CCS with joint recovery on Azure batch. This software is based on [JUDI](https://github.com/slimgroup/JUDI.jl), which uses the highly optimized time-domain finite-difference propagators of [Devito](https://github.com/devitocodes/devito) and explores the serverless parallelism through [JUDI4Cloud](https://github.com/slimgroup/JUDI4Cloud.jl).

## TO-DO

First, please put your Azure credential file in this directory (`Seis4CCS.jl/credentials.json`). e.g. do

```bash
touch credentials.json
vi credentials.json
```

and copy-paste your credentials to this file.

Then, under `Seis4CCS.jl`, run

```bash
L=a nsrc=b vm=c nth=d niter=e bs=f snr=g gamma=h bash runall.sh
```

where integer `a` is the number of vintages in the experiment (maximum is 879), integer `b` is number of jittered sample seismic sources in each vintage (different for each vintage, suggested to be more than `8`), `c` is the type of virtual machine for seismic data simulation and imaging, `d` is the number of threads on each computing node in Azure batch, `e` is the number of joint recovery iterations, `f` is the batchsize in a JRM iteration, float `g` is the SNR of the noise to be added in the data, float `h` is the weight `gamma` on common component (suggested to be `a/2`). The default values of `a,b,c,d,e,f,g,h` are `2,8,Standard_E8s_v3,4,16,4,0,1.0`.

When the aformentioned program finishes, the generated seismic data will be saved in `./data` folder and each iteration of joint recovery results will be saved in `./results` folder. To test if everything runs without error in a small scale, you could do

```bash
bash runall.sh
```

These scripts are suggested to run from docker image, which will be provided upon request.

Notice that if you set `L`, `nsrc`, `niter` to be very large, it may require massive disk space and memory on your local machine. It is also suggested to add `--oom-kill-disable` when you run the docker image.

There is `rclone` installed in the docker image. Feel free to use it to upload the results to your preferred online platforms.

## Authors

If any question, please contact

Ziyi (Francis) Yin, ziyi.yin@gatech.edu    
Mathias Louboutin, mlouboutin3@gatech.edu    
Felix J. Herrmann, felix.herrmann@gatech.edu     

This research was carried out with the support of Georgia Research Alliance and partners of the ML4Seismic Center.