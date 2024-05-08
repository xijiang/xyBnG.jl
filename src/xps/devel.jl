using Mmap
using Serialization
using xyBnG.XY
using xyBnG.xyTypes: Trait
using xyBnG.Breeding: phenotype!, predict!

function devel(bar::AbstractString)
    @info "Directional selection test"
    milk   = Trait("milk",   0.3, 10000; sex = 0, vd = 0.1)
    growth = Trait("growth", 0.5, 10000; sex = 2, vd = 0.2)
    dir = "rst/$bar"
    lmp = deserialize("$dir/$bar.lmp")
    ped = deserialize("$dir/$bar.ped")
    predict!(1:nrow(ped), ped, milk, growth)
end
