module Founder

using DataFrames
using Distributions
using LinearAlgebra
using Mmap
using Random
using Serialization
using xyBnG.Conn.TS
using xyBnG.XY
import xyBnG.xyTypes: Cattle, Species, Trait

include("fdr/xybase.jl")
include("fdr/extern.jl")
include("fdr/sample.jl")

end # module Founder
