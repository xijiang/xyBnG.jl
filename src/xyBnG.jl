module xyBnG

include("misc/types.jl")
include("misc/utils.jl")
include("xy/xy.jl")
include("conn/conn.jl")
include("founder/founder.jl")
include("breeding/breeding.jl")
include("relationships/relationship.jl")
include("summary/summary.jl")
include("xps/xps.jl")

import xyBnG.xyTypes: Species, Cattle, Trait, Plan
import xyBnG.Founder: ts_base

export Species, Cattle, Trait, Plan, ts_base

end # module xyBnG
