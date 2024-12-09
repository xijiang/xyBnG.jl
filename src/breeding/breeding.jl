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
using Statistics
import StatsBase: sample, Weights
using xyBnG.XY
using xyBnG.xyTypes: Trait, Cattle, Species, Plan, name

include("utils.jl")
include("phenotype.jl")
include("predict.jl")
include("select.jl")
include("reproduce.jl")
include("ocs.jl")
include("xecodes.jl")

end # module Breeding
