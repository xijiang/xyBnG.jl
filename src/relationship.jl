module RS

using DataFrames
using LinearAlgebra
using Mmap
using Statistics
using xyBnG.XY
import xyBnG.Util: memavail

include("rs/grm.jl")
include("rs/irm.jl")
include("rs/nrm.jl")
include("rs/phasedibd.jl")

end # module RS
