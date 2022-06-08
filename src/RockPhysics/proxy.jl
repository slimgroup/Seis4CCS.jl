export VtoK, Ktoϕ

function VtoK(v::AbstractMatrix{T}, d::Tuple{T, T}; α=1.03f0) where T

    n = size(v)
    idx_wb = find_water_bottom(v.-minimum(v))
    idx_ucfmt = find_water_bottom((v.-3.5f0).*(v.>3.5f0))
    Kh = zeros(Float32, n)
    capgrid = Int(round(50f0/d[2]))
    for i = 1:n[1]
        Kh[i,1:idx_wb[i]-1] .= 1f-10  # water layer
        Kh[i,idx_wb[i]:idx_ucfmt[i]-capgrid-1] .= α*v[i,idx_wb[i]:idx_ucfmt[i]-capgrid-1].+15f0  # secondary seal
        Kh[i,idx_ucfmt[i]-capgrid:idx_ucfmt[i]-1] .= 1f-3   # seal
        Kh[i,idx_ucfmt[i]:end] .= α*v[i,idx_ucfmt[i]:end].+200f0    # reservoir
    end
    return Kh
end

function Ktoϕ(K::AbstractMatrix{T}, d::Tuple{T, T}) where T

    n = size(K)
    idx_ucfmt = find_water_bottom((K.-1f-3).*((K.==1f-3).|(K.>200f0)))
    capgrid = Int(round(50f0/d[2])) # 50m seal
    phi = zeros(Float32,n)
    for i = 1:n[1]
        for j = 1:n[2]
            p = Polynomial([-0.0314^2*K[i,j],2*0.0314^2*K[i,j],-0.0314^2*K[i,j],1.527^2])
            phi[i,j] = minimum(real(roots(p)[findall(real(roots(p)).== roots(p))]))
        end
        phi[i,idx_ucfmt[i]-capgrid:idx_ucfmt[i]-1] = Float32.(range(0.056,stop=0.1,length=capgrid)) # porosity gradually decreases at the seal
    end
    return phi
end