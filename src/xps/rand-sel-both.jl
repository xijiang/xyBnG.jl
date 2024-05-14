"""
    rand_sel_both()

Ramdom seleciton on trait milk, and growth. The trait's EBV were randomly
assigned. The weight are 0.6 and 0.4.
"""
function rand_sel_both()
    @info "Directional selection test: index selection"
    foo = "cattle"
    dir, bar = "rst/$foo", "test"
    # backup the founder files
    for ext in ["lmp", "ped", "xy"]
        cp("$dir/$foo.$ext", "$dir/$bar.$ext", force = true)
    end
    milk   = Trait("milk",   0.3, 10000; sex = 0, vd = 0.1)
    growth = Trait("growth", 0.5, 10000; sex = 2, vd = 0.2)

    lmp = deserialize("$dir/$bar.lmp")
    ped = deserialize("$dir/$bar.ped")
    xy  = "$dir/$bar.xy"
    plan = Plan(20, 40, 200) # difault :hierarchical
    wgt = Dict{String, Float64}("milk" => 0.6, "growth" => 0.4)

    for igrt in 1:5
        println()
        @info lpad("<--- Generation $igrt/5 --->", 40)
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, milk, growth)
        predict!(ids, ped, milk, growth)
        ng = Select(ids, plan, ped, wgt)
        reproduce!(ng, ped, xy, lmp, milk, growth)
    end
    serialize("$dir/$bar.ped", ped)
end