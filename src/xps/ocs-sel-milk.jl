function ocs_sel_milk(; bar = "ocs")
    dir, fixed = "rst/cattle", ["grt"]
    isfile("$dir/rand.xy") || rand_sel_milk(bar = bar)

    for ext in ["ped", "xy"]
        cp("$dir/rand.$ext", "$dir/$bar.$ext")
    end

    @info "Directional selection test: index selection"
    milk   = Trait("milk",   0.3, 10000; sex = 0, vd = 0.1)
    growth = Trait("growth", 0.5, 10000; sex = 2, vd = 0.2)
    lmp = deserialize("$dir/$bar.lmp")
    ped = deserialize("$dir/$bar.ped")
    xy  = "$dir/$bar.xy"
    plan = Plan(20, 40, 200) # difault :hierarchical

    for igrt in 1:1 # AABLUP
        println()
        @info lpad("<--- Generation $igrt/10 --->", 40)
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, milk, growth)
        giv, R = begin
            G = nrm(ped)
            giv = inv(G)
            R = G[ids, ids]
            giv, R
        end
        predict!(ids, ped, fixed, giv, milk)
        tgt = select(ped, "ebv_milk" => "ebv", "sex")[ids, :]
        K = 0.038 # ≈ 1/Ne
        c = TM1997(tgt, R, K)
        return c
        #ng = Select(ids, c, ped, milk)
        #reproduce!(ng, ped, xy, lmp, milk, growth)
    end
    serialize("$dir/$bar.ped", ped)
end


## ToDo
# 1. unify Select
# 2. finish this function
# 3. Dos ocs.
# 4. More simulation on dos ocs
# 5. Oda's scheme.