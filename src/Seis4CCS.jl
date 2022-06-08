# Module with functions for flow modeling and rock physics
# Author: Ziyi Yin, ziyi.yin@gatech.edu
# Date: June, 2022

__precompile__()

module Seis4CCS

export Seis4CCSPATH
Seis4CCSPATH = dirname(pathof(Seis4CCS))

include("FlowSimulation/FlowSimulation.jl")
include("RockPhysics/RockPhysics.jl")


end