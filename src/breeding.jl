"""
    module Breeding

## Provides
- `function sumMap(lmp; cM = 1e6)`
- `function crossovers(lms)`
- `function gamete(prt, hap, lms)`
"""
module Breeding

using DataFrames
using Distributions
using Mmap
using Random
using Serialization
using xyBnG.XY
using xyBnG.xyTypes

include("brd/utils.jl")
include("brd/phenotype.jl")
include("brd/predict.jl")
include("brd/select.jl")
include("brd/reproduce.jl")
include("brd/ocs.jl")

end # module Breeding