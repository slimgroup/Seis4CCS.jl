export FlowSimulation

module FlowSimulation

    using FwiFlow
    @warn "This module will be replaced by a straightforward non-tensorflow implementation that scales to 3D in the future, e.g. based on the OPM/GEOSX"
    include("FlowSimulationFunctions.jl") 

end