export Patchy, find_water_bottom

#### Patchy saturation model adapted from FwiFlow.jl

function Patchy(sw::AbstractMatrix{T}, vp::AbstractMatrix{T}, rho::AbstractMatrix{T}, phi::AbstractMatrix{T}, d::Tuple{T, T}; bulk_fl1 = 2.735f9, bulk_fl2 = 0.125f9,ρw = 7.766f2, ρo = 1.053f3) where T

    n = size(sw)
    capgrid = Int(round(50f0/d[2]))
    vp = vp * 1f3
    vs = vp ./ sqrt(3f0)
    idx_wb = maximum(find_water_bottom(vp.-minimum(vp)))
    idx_ucfmt = find_water_bottom((vp.-3500f0).*(vp.>3500f0))

    bulk_sat1 = rho .* (vp.^2f0 - 4f0/3f0 .* vs.^2f0) * 1f3
    shear_sat1 = rho .* (vs.^2f0) * 1f3

    bulk_min = zeros(Float32,size(bulk_sat1))

    bulk_min[findall(vp.>=3500f0)] .= 5f10   # mineral bulk moduli
    bulk_min[findall(vp.<3500f0)] .= 1.2f0 * bulk_sat1[findall(vp.<3500f0)] # mineral bulk moduli

    patch_temp = bulk_sat1 ./(bulk_min .- bulk_sat1) - bulk_fl1 ./ phi ./ (bulk_min .- bulk_fl1) + bulk_fl2 ./ phi ./ (bulk_min .- bulk_fl2)

    for i = 1:n[1]
        patch_temp[i,idx_ucfmt[i]-capgrid:idx_ucfmt[i]-1] = patch_temp[i,idx_ucfmt[i]:idx_ucfmt[i]+capgrid-1]
    end

    bulk_sat2 = bulk_min./(1f0./patch_temp .+ 1f0)
    bulk_sat2[findall(bulk_sat2-bulk_sat1.>0)] = bulk_sat1[findall(bulk_sat2-bulk_sat1.>0)]

    bulk_new = (bulk_sat1+4f0/3f0*shear_sat1).*(bulk_sat2+4f0/3f0*shear_sat1) ./( (1f0.-sw).*(bulk_sat2+4f0/3f0*shear_sat1) 
    + sw.*(bulk_sat1+4f0/3f0*shear_sat1) ) - 4f0/3f0*shear_sat1

	bulk_new[:,1:idx_wb] = bulk_sat1[:,1:idx_wb]
    bulk_new[findall(sw.==0)] = bulk_sat1[findall(sw.==0)]
    rho_new = rho + phi .* sw * (ρw - ρo) / 1f3
    rho_new[findall(sw.==0)] = rho[findall(sw.==0)]
    Vp_new = sqrt.((bulk_new+4f0/3f0*shear_sat1)./rho_new/1f3)
    Vp_new[findall(sw.==0)] = vp[findall(sw.==0)]

    return Vp_new/1f3, rho_new
end
function Patchy(sw::AbstractArray{T, 3}, vp::AbstractMatrix{T}, rho::AbstractMatrix{T}, phi::AbstractMatrix{T}, d::Tuple{T, T}; bulk_fl1 = 2.735f9, bulk_fl2 = 0.125f9,ρw = 7.766f2, ρo = 1.053f3) where T

    stack = [Patchy(sw[i,:,:], vp, rho, phi, d; bulk_fl1=bulk_fl1, bulk_fl2=bulk_fl2, ρw = ρw, ρo=ρo) for i = 1:size(sw,1)]
    return [stack[i][1] for i = 1:size(sw,1)], [stack[i][2] for i = 1:size(sw,1)]
end

#### verbatim copy from JUDI.jl

function find_water_bottom(m::AbstractArray{avDT,2};eps = 1e-4) where {avDT}
    #return the indices of the water bottom of a seismic image
    n = size(m)
    idx = zeros(Integer, n[1])
    for j=1:n[1]
        k=1
        while true
            if abs(m[j,k]) > eps
                idx[j] = k
                break
            end
            k += 1
        end
    end
    return idx
end

find_water_bottom(m::AbstractArray{avDT,3};eps = 1e-4) where {avDT} = reshape(find_water_bottom(reshape(m, :, size(m)[end]);eps=eps), size(m,1), size(m,2))
