module xyBnG

include("types.jl")
include("utils.jl")
include("xy.jl")
include("conn.jl")
include("founder.jl")
include("breeding.jl")
include("relationship.jl")
include("summary.jl")
include("xps.jl")

import xyBnG.xyTypes: Species, Cattle, Trait, Plan
import xyBnG.Founder: ts_base

export Species, Cattle, Trait, Plan, ts_base

end # module xyBnG
