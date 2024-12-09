"""
    Sum
The ``Sum`` module provides functions to calculate summary statistics of the
simulations. It has majorly two parts. One is to summarize each repeat of a
simulation. The other contains a few notebooks for the presentation of the
summaries.
"""
module Sum

using DataFrames
using LinearAlgebra
using Serialization
using Statistics
using xyBnG.XY
using xyBnG.RS
import xyBnG.xyTypes: Cattle, Species, Trait, Plan, name

include("savepar.jl")
include("tibd.jl")
include("ffcv.jl")
include("xysum.jl")
include("savesum.jl")
include("cormat.jl")
include("snphet.jl")
include("fgrm.jl")
include("resum.jl")

end # module Sum
