module Xps

using DataFrames
using Mmap
using Serialization
import xyBnG.xyTypes: Cattle, Species, Trait
import xyBnG.Founder: sim_base, sample_founder
import xyBnG.Conn.TS: sample2xy

include("xps/samples.jl")

end # module Xps