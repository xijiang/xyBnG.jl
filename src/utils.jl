"""
    module Util

Some utility functions.

## Provides
- `function memavail()`
- `function blksz(n::Int, m::Int)`
- `function commas(num::T) where T <: Integer`
- `function Ceiling(; nQTL = 1000, nRepeat = 1000, ad = Normal(), qd =
  Beta(0.75, 0.75))`
"""
module Util
using Distributions

"""
    memavail()

Return the size of available memory on a Linux alike system. Or, the free memory
size is returned.
"""
function memavail()
    if Sys.isunix()
        for line in eachline("/proc/meminfo")
            if occursin("MemAv", line)
                return parse(Int, split(line)[2]) * 1024
            end
        end
    else
        return Sys.free_memory()
    end
end

"""
    blksz(n::Int, m::Int)

Group `n` ID as evenly as possible with each group has maximally `m` ID. The
function returns the group size.
"""
function blksz(n::Int, m::Int)
    n % m == 0 ? m : begin
        ng = Int(ceil(n / m))
        n % ng == 0 ? n ÷ ng : n ÷ ng + 1
    end
end

"""
    commas(num::T)

Display `num` with commas.
"""
function commas(num)
    str = string(num)
    return replace(str, r"(?<=[0-9])(?=(?:[0-9]{3})+(?![0-9]))" => ",")
end

"""
    Ceiling(; nQTL = 1000, nRepeat = 1000, ad = Normal(), qd = Beta(0.75, 0.75))

Calculate the ceiling height of a population with `nQTL` QTLs, whose effects are
drawn from `ad` and allele frequencies are drawn from `qd`.
"""
function Ceiling(;
    nQTL = 1000,
    nRepeat = 1000,
    ad = Normal(),         # can be of other distributions
    qd = Beta(0.75, 0.75), # allele frequency distribution
)
    C = 0.0        # Ceiling height
    for _ in 1:nRepeat
        a = rand(ad, nQTL) # effects of SNP 1s
        p = rand(qd, nQTL) # U-shaped allele frequencies
        σ = sqrt(sum(2p .* (1 .- p) .* a .^ 2)) # genic std
        C += 2sum(a[a.>0.0]) / σ
    end
    C /= nRepeat
end

end # module util
