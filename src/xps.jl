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
using xyBnG.Sum
import xyBnG.xyTypes: Cattle, Species, Trait, Plan, name
import xyBnG.Breeding: phenotype!, Predict!, Select, reproduce!, TM1997, TM2024
import xyBnG.Founder: sim_base, sample_founder, uniq
import xyBnG.Conn.TS: sample2xy
import xyBnG.RS: nrm, irm, grm, xirm

include("xps/a-tale-of-selction.jl")
include("xps/one-trait.jl")
include("xps/devel.jl")
include("xps/dosblup.jl")

end # module xps