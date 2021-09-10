function ContJitter(n, d, num)

    base_obn = convert(Array{Float32},range(d[1],stop=(n[1]-1)*d[1],length = num)) # spacing = 250m on average
    avg_spc = (n[1]-2)*d[1]/(num-1)
    OBN_loc = zeros(Float32,num)

    for i = 1:num
        OBN_loc[i] = max(min(base_obn[i] + avg_spc*rand(Float32) - avg_spc/2, (n[1]-1)*d[1]),0f0)
    end

    return OBN_loc
end

function Jitt_function(N, gamma, xi)
    # 2020.12.09 Yijun
    # jittered subsampling function for simple test 
    # N: the full size 
    # gamma: undersampling factor
    # xi: control the maximum gap between the neighbor subsampling indexes
    # gamma=4 --> 75 % subsampling
    #rng(rseed);
    test1 = rand(1,Int64(round(N/gamma))).*(xi-1) .- (xi-1)/2 ;
    test2 = ones(1, Int64(round(N/gamma)));
    test3 = cumsum(test2; dims=2); 
    n = Int64.(round.((1-gamma)/2 .+ gamma .* test3 .+ test1));
    for i = 1:length(n)
        if n[i] > N
            n[i] = N
        end
    end
    return n
end