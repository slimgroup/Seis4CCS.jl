# Module with functions for flow modeling and rock physics
# Author: Ziyi Yin, ziyi.yin@gatech.edu
# Date: June, 2022

__precompile__()

module Seis4CCS

using JOLI, JUDI, FwiFlow

export Seis4CCSPATH
Seis4CCSPATH = dirname(pathof(Seis4CCS))

# submodule Simulation
include("Simulation/Simulation.jl")

# submodule Imaging
include("Imaging/Imaging.jl")


end