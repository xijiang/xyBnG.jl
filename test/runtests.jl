using xyBnG
using Test
using DataFrames

import CompatHelperLocal as CHL
CHL.@check()

include("tst-xy.jl")
include("tst-founder.jl")
include("tst-rs.jl")
