export flow, comp_grid

const K_CONST =  9.869232667160130e-16 * 86400 * 1e3
const ALPHA = 1.0
const SRC_CONST = 86400.0 #
const GRAV_CONST = 9.8    # gravity constant

## Only 2D flow simulation is supported now
mutable struct comp_grid       # set up the grid structure
    n  :: Tuple{Integer, Integer} # x, z
    h  :: Float64 # meter
    hy :: Float64 # meter
    nt :: Integer 
    dt :: Float64  # day
end

function comp_grid(n::Tuple{Integer, Integer}, h::Tuple{T1, T1}, hy::T2, nt::Integer, dt::T3) where {T1, T2, T3}
    @assert h[1] == h[2] "Now we only support equal spacing in x and z direction"
    return comp_grid(n, Float64(h[1]), Float64(hy), nt, Float64(dt))
end

function comp_grid(n::Tuple{Integer, Integer}, h::T1, hy::T2, nt::Integer, dt::T3) where {T1, T2, T3}
    return comp_grid(n, h[1], Float64(hy), nt, Float64(dt))
end

## Simulation structure for FwiFlow
mutable struct Ctx
  m; n; h; NT; Δt; Z; X; ρw; ρo;
  μw; μo; K; g; ϕ; qw; qo; sw0
end

function tfCtxGen(m,n,h,NT,Δt,Z,X,ρw,ρo,μw,μo,K,g,ϕ,qw,qo,sw0)
    tf_h = constant(h)
    tf_Δt = constant(Δt)
    tf_Z = constant(Z)
    tf_X= constant(X)
    tf_ρw = constant(ρw)
    tf_ρo = constant(ρo)
    tf_μw = constant(μw)
    tf_μo = constant(μo)
    tf_K = constant(K)
    tf_g = constant(g)
    tf_ϕ = constant(ϕ)
    tf_qw = constant(qw)
    tf_qo = constant(qo)
    tf_sw0 = constant(sw0)
    return Ctx(m,n,tf_h,NT,tf_Δt,tf_Z,tf_X,tf_ρw,tf_ρo,tf_μw,tf_μo,tf_K,tf_g,tf_ϕ,tf_qw,tf_qo,tf_sw0)
end

# variables : sw, u, v, p
# (time dependent) parameters: qw, qo, ϕ
function onestep(sw, p, m, n, h, Δt, Z, ρw, ρo, μw, μo, K, g, ϕ, qw, qo)
    # step 1: update p
    λw = sw.*sw/μw
    λo = (1-sw).*(1-sw)/μo
    λ = λw + λo
    q = qw + qo + λw/(λo+1e-16).*qo
    potential_c = (ρw - ρo)*g .* Z

    # Step 1: implicit potential
    Θ = upwlap_op(K * K_CONST, λo, potential_c, h, constant(0.0))

    load_normal = (Θ+q/ALPHA) - mean(Θ+q/ALPHA)

    p = upwps_op(K * K_CONST, λ, load_normal, p, h, constant(0.0), constant(0)) # potential p = pw - ρw*g*h 

    # step 2: implicit transport
    sw = sat_op(sw, p, K * K_CONST, ϕ, qw, qo, μw, μo, sw, Δt, h)
    return sw, p
end


function imseq(tf_ctx)
    ta_sw, ta_p = TensorArray(tf_ctx.NT+1), TensorArray(tf_ctx.NT+1)
    ta_sw = write(ta_sw, 1, tf_ctx.sw0)
    ta_p = write(ta_p, 1, constant(zeros(tf_ctx.m, tf_ctx.n)))
    i = constant(1, dtype=Int32)
    function condition(i, tas...)
        i <= tf_ctx.NT
    end
    function body(i, tas...)
        ta_sw, ta_p = tas
        sw, p = onestep(read(ta_sw, i), read(ta_p, i), tf_ctx.m, tf_ctx.n, tf_ctx.h, tf_ctx.Δt, tf_ctx.Z, tf_ctx.ρw, tf_ctx.ρo, tf_ctx.μw, tf_ctx.μo, tf_ctx.K, tf_ctx.g, tf_ctx.ϕ, tf_ctx.qw[i], tf_ctx.qo[i])
        ta_sw = write(ta_sw, i+1, sw)
        ta_p = write(ta_p, i+1, p)
        i+1, ta_sw, ta_p
    end

    _, ta_sw, ta_p = while_loop(condition, body, [i, ta_sw, ta_p])
    out_sw, out_p = stack(ta_sw), stack(ta_p)
end

function flow(K, ϕ, qw_value, qo_value, grid;
    ρw=501.9, ρo=1053.0, μw=0.1, μo=1.0, g=GRAV_CONST, S0=nothing)
    qw = tf.placeholder(tf.float64)
    qo = tf.placeholder(tf.float64)
    if isnothing(S0)
        S0 = zeros(grid.n[2], grid.n[1])
    end
    qw_value = permutedims(qw_value, [1, 3, 2]) * (1/grid.h^2)/grid.hy * SRC_CONST
    qo_value = permutedims(qo_value, [1, 3, 2]) * (1/grid.h^2)/grid.hy * SRC_CONST
    n = size(K)
    K = K'
    ϕ = ϕ'
    X = reshape(repeat((1:n[1])*grid.h, outer = n[2]), n[1], n[2])'
    Z = reshape(repeat((1:n[2])*grid.h, outer = n[1]), n[2], n[1])
    tfCtxTrue = tfCtxGen(grid.n[2],grid.n[1],grid.h,grid.nt,grid.dt,Z,X,ρw,ρo,μw,μo,K,g,ϕ,qw,qo,S0)
    sp = imseq(tfCtxTrue)
    sess = Session(); init(sess)
    S, p = sess.run(sp, Dict(qw=>qw_value,qo=>qo_value))
    return permutedims(S, [1, 3, 2]), permutedims(p, [1, 3, 2])
end 