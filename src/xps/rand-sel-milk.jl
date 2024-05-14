"""
    rand_sel_milk()

Ramdom seleciton on trait milk. The trait's EBV are randomly assigned.
"""
function rand_sel_milk(; bar = "rand")
    @info "Directional selection test"
    foo, dir = "base", "rst/cattle"
    # backup the founder files
    for ext in ["ped", "xy"]
        cp("$dir/$foo.$ext", "$dir/$bar.$ext", force = true)
    end
    milk   = Trait("milk",   0.3, 10000; sex = 0, vd = 0.1)
    growth = Trait("growth", 0.5, 10000; sex = 2, vd = 0.2)

    lmp = deserialize("$dir/$foo.lmp")
    ped = deserialize("$dir/$bar.ped")
    xy  = "$dir/$bar.xy"
    plan = Plan(20, 40, 200) # default :hierarchical

    for igrt in 1:5
        println()
        @info lpad("<--- Generation $igrt/5 --->", 40)
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, milk, growth)
        predict!(ids, ped, milk, growth)
        ng = Select(ids, plan, ped, milk)
        reproduce!(ng, ped, xy, lmp, milk, growth)
    end
    serialize("$dir/$bar.ped", ped)
end
