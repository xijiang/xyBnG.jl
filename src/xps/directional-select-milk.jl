"""
    directional_select_milk()

Directional selection test. This selection is on EBV of trait `milk`. The ID in
the last generation are as candidates. The example here is ABLUP. If we replace
`nrm` here with `grm`, we will have a GBLUP example. Note the argument of
`fixed` effects should be a vector of strings. These strings are the names in
the pedigree `ped`.
"""
function directional_select_milk(; bar = "ablup")
    dir, fixed = "rst/cattle", ["grt"]
    isfile("$dir/rand.xy") || rand_sel_milk(bar = "rand")

    for ext in ["ped", "xy"]
        cp("$dir/rand.$ext", "$dir/$bar.$ext")
    end

    @info "Directional selection test: index selection"
    milk   = Trait("milk",   0.3, 10000; sex = 0, vd = 0.1)
    growth = Trait("growth", 0.5, 10000; sex = 2, vd = 0.2)
    lmp = deserialize("$dir/base.lmp")
    ped = deserialize("$dir/$bar.ped")
    xy  = "$dir/$bar.xy"
    plan = Plan(20, 40, 200) # default :hierarchical

    for igrt in 1:10
        println()
        @info lpad("<--- Generation $igrt/10 --->", 40)
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, milk, growth)
        G = nrm(ped)
        giv = inv(G)
        predict!(ids, ped, fixed, giv, milk)
        ng = Select(ids, plan, ped, milk)
        reproduce!(ng, ped, xy, lmp, milk, growth)
    end
    serialize("$dir/$bar.ped", ped)
end
