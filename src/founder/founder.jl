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

include("fwdsim.jl")
include("xybase.jl")
include("extern.jl")
include("sample.jl")

end # module Founder
