"""
    module xps

Some examples to use package `xyBnG`. Note these functions depend on the results
of previous functions. Hence, to run these functions, one must have run the
functions before them.
"""
module xps

using DataFrames
using LinearAlgebra
using Markdown
using Mmap
using Serialization
using Term
using xyBnG.XY
import xyBnG.xyTypes: Cattle, Species, Trait, Plan, name
import xyBnG.Breeding: phenotype!, predict!, Select, reproduce!, TM1997
import xyBnG.Founder: sim_base, sample_founder
import xyBnG.Conn.TS: sample2xy
import xyBnG.RS: nrm

include("xps/directional-select-milk.jl")
include("xps/ocs-sel-milk.jl")
include("xps/a-tale-of-selction.jl")

end # module xps