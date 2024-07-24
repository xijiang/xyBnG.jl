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
using LinearAlgebra
using Mmap
using Random
using Serialization
using SparseArrays
import StatsBase: sample, Weights
using xyBnG.XY
using xyBnG.xyTypes: Trait, Cattle, Species, Plan, name

include("brd/utils.jl")
include("brd/phenotype.jl")
include("brd/predict.jl")
include("brd/select.jl")
include("brd/reproduce.jl")
include("brd/ocs.jl")
include("brd/xecodes.jl")

end # module Breeding
