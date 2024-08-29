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
using Random
using Serialization
using Term
using xyBnG.XY
using xyBnG.Conn
using xyBnG.Sum
import xyBnG.xyTypes: Cattle, Species, Trait, Plan, name
import xyBnG.Breeding: phenotype!, Predict!, Select, reproduce!, TM1997, TM2024
import xyBnG.Founder: ts_base, uniq, sample_xy, macs_base
import xyBnG.RS: nrm, irm, grm, xirm

# Schemes
include("xps/schemes/cull.jl")
include("xps/schemes/ocs-2sex.jl")
include("xps/schemes/ocs-male.jl")

# Examples
include("xps/chkbase.jl")
include("xps/init.jl")
include("xps/a-tale-of-selction.jl")
include("xps/one-trait.jl")
include("xps/dosblup.jl")
include("xps/o3cs.jl")
include("xps/oc-sires.jl")
include("xps/mini.jl")

end # module xps
