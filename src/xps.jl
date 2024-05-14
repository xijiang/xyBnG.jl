"""
    module xps

Some examples to use package `xyBnG`. Note these functions depend on the results
of previous functions. Hence, to run these functions, one must have run the
functions before them.
"""
module xps

using DataFrames
using LinearAlgebra
using Mmap
using Serialization
using xyBnG.XY
import xyBnG.xyTypes: Cattle, Species, Trait, Plan, name
import xyBnG.Breeding: phenotype!, predict!, Select, reproduce!
import xyBnG.Founder: sim_base, sample_founder
import xyBnG.Conn.TS: sample2xy
import xyBnG.RS: nrm

include("xps/sim-cattle-base.jl")
include("xps/sample-a-founder.jl")
include("xps/rand-sel-milk.jl")
include("xps/rand-sel-both.jl")
include("xps/directional-select-milk.jl")

end # module xps