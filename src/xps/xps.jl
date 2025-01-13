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
using Statistics
using xyBnG.XY
using xyBnG.Conn
using xyBnG.Sum
import xyBnG.xyTypes: Cattle, Species, Trait, Plan, name
import xyBnG.Breeding: phenotype!, Predict!, Select, reproduce!, TM1997, TM2024
import xyBnG.Founder: ts_base, uniq, sample_xy, macs_base
import xyBnG.RS: nrm, irm, grm, xirm

# Schemes
include("schemes/cull.jl")
include("schemes/ocs-2sex.jl")
include("schemes/ocs-male.jl")
include("schemes/file-irm.jl")

# Examples
include("chkbase.jl")
include("init.jl")
include("a-tale-of-selction.jl")
include("one-trait.jl")
include("dosblup.jl")
include("o3cs.jl")
include("oc-sires.jl")
include("mini.jl")
include("test-scheme.jl")
include("f0.jl")

end # module xps
