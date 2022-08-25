# Seis4CCS

[Seis4CCS] is a framework for [simulation-based seismic monitoring design](https://slim.gatech.edu/research/geological-carbon-storage#simulation-based-monitoring-design) for [geological carbon storage](https://slim.gatech.edu/research/geological-carbon-storage) (GCS) monitoring. At [SLIM], we aim to reduce the seismic operating costs by optimizing acquisition design and to answer questions such as

- How often and when (early and/or late in a GCS projects) do seismic surveys need to be performed to mitigate long-term risks?
- How can the sensitivity of seismic monitoring systems be improved while keeping costs down?
- To what degree do time-lapse surveys have to be replicated to achieve high degrees of repeatability?

To answer these questions and to help drive innovations in seismic monitoring acquisition design and imaging, we are developing an open-source software platform [Seis4CCS] that allows users to perform simulation-based monitoring design and test novel time-lapse acquisition and imaging technologies *in silico* at scale.

# Installation

Seis4CCS can be installed with the standard julia package manager:

```julia
] add https://github.com/slimgroup/Seis4CCS.jl.git
```

or in a developer mode

```julia
] dev https://github.com/slimgroup/Seis4CCS.jl.git
```

# Modules

This repository currently contains two modules, namely

- flow simulation, which builds a wrapper for two-phase flow simulations based on [FwiFlow]
- rock physics, which converts time-varying CO~2~ concentration to wave properties (e.g. wavespeed and density) based on [patchy saturation model](https://www.cambridge.org/core/books/quantitative-seismic-interpretation/EB6A36B78CCF07187723F6F5364EDCF8)

The third module in the simulation-based monitoring is wave physics and it can be simulated in [JUDI], which uses the highly optimized time-domain finite-difference propagators of [Devito]. Instruction for JUDI installation can be found in JUDI's README.

# Notebooks

We provide three notebooks for flow simulation, rock physics and seismic imaging respectively in the `notebooks` folder.

# Author

Ziyi Yin, ziyi.yin@gatech.edu

[Seis4CCS]:https://github.com/slimgroup/Seis4CCS.jl
[SLIM]:https://slim.gatech.edu/
[FwiFlow]:https://github.com/lidongzh/FwiFlow.jl
[JUDI]:https://github.com/slimgroup/JUDI.jl
[Devito]:https://www.devitoproject.org/