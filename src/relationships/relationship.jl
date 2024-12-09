module RS

using DataFrames
using LinearAlgebra
using Mmap
using Statistics
using xyBnG.XY
import xyBnG.Util: memavail

include("grm.jl")
include("irm.jl")
include("nrm.jl")
include("phasedibd.jl")

end # module RS
