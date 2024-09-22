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
import xyBnG.xyTypes: Trait

include("sum/savepar.jl")
include("sum/tibd.jl")
include("sum/ffcv.jl")
include("sum/xysum.jl")
include("sum/savesum.jl")
include("sum/cormat.jl")
include("sum/snphet.jl")
include("sum/fgrm.jl")
include("sum/resum.jl")

end # module Sum
