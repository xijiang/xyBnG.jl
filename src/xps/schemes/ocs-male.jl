
function aaom(test, foo, bar, lmp, ngn, trait, fixed, plan, dF, F0)
    @info " - Directional selection FXOCS for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force = true)

    for ign = 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        G = nrm(ped)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        g22 = G[ids, ids]
        ng = Select(ids, plan, ped, g22, trait, dF, ign; F0 = F0)
        reproduce!(ng, ped, xy, lmp, trait)
    end
end
