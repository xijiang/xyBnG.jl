module xyBnG

include("types.jl")
include("utils.jl")
include("xy.jl")
include("conn.jl")
include("founder.jl")
include("breeding.jl")
include("relationship.jl")
include("xps.jl")

import xyBnG.xyTypes: Species, Cattle, Trait
import xyBnG.Founder: sim_base
import xyBnG.xps: devel

export Species, Cattle, Trait, sim_base
export devel

end # module xyBnG
