"""
    fwdsim(nn::Matrix{Int}; chr::Float64 = 1.0, mr::Float64 = 3.0, bppm::Int = 100_000_000)
This is a reproduction of my `pgsnp.cpp` written in 2010. This Julia version
removes the constrains about stack size. It also added variation of Ne in
history.
"""
function fwdsim(
    nn::Matrix{Int};    # Ne and their generations
    chr::Float64 = 1.0, # chromosome length in Morgans
    mr::Float64 = 3.0,  # mutation rate per Morgan per meiosis
    bppm::Int = 100_000_000, # base pairs per Morgan
)
    @info "Forward simulation of:
          chromosome length $chr Morgans,
          mutation rate $mr, 
          $bppm base pairs per Morgan,
          in a history ([ne, n-grt;...])
          $nn
          "

    # Check parameters
    @assert 1e-5 < chr < 6.0 "chromosome length out of range (1e5, 6.0)"
    @assert 1e-2 < mr < 7.0 "mutation rate out of range (1e-2, 7.0)"
    @assert 1e5 < bppm < 1e9 "base pairs per Morgan out of range (1e5, 1e9)"
    @assert size(nn, 2) == 2 "nn must be a matrix with 2 columns"
    @assert all(nn[:, 1] .> 4) "ne must be greater than 4"
    @assert all(nn[:, 2] .> 2) "n-grt must be greater than 0"

    mm = 120_000   # maximum number of mutations
    mtv = 240_000
    cp = 1_000     # period to comb out fixed loci
    tp = 500       # progress period
    tvec = zeros(Int, mtv)

    parents = zeros(Int, nn[1, 1], 2, mm)
    children = zeros(Int, nn[1, 1], 2, mm)
    nmp = zeros(Int, nn[1, 1], 2)
    nmc = zeros(Int, nn[1, 1], 2)
    λₘ = chr * mr
    #mx = mr < 1 ? 

    for (ne, ngrt) in eachrow(nn)
        @info "Simulating $ne individuals for $ngrt generations"
    end
end
