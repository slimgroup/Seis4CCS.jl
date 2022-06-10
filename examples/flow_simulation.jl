using Seis4CCS.FlowSimulation

### This script follows the tutorial at https://lidongzh.github.io/FwiFlow.jl/dev/tutorials/flow/

using SlimPlotting, PyPlot

# Hyperparameter for flow simulation
n = (30, 15)    # num of cells
h = 30.0        # meter
nt = 50       # number of time steps
dt = 20.0       # day

grid_ = comp_grid(n, h, h, nt, dt);

K = 20.0 .* ones(n) # permeability 
ϕ = 0.25 .* ones(n) # porosity

qw = zeros(nt, n[1], n[2]);
qw[:,3,9] .= 0.005;   # in [m^3/s]
qo = zeros(nt, n[1], n[2]);
qo[:,28,9] .= -0.005;

@time S, p = flow(K, ϕ, qw, qo, grid_);

figure(figsize=(20,12))
for i = 1:6
    subplot(2,3,i)
    plot_velocity(S[10*i-9,:,:]', (h, h); new_fig=false, vmax=1)
    colorbar()
    title("snapshot at day $((i-1)*10*dt)")
end