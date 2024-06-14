function devel(; op = 1)
    if op == 1
        one_trait()
    elseif op == 2
        tstDir = "rst/onetrt"
        milk = Trait("milk", 0.25, 10000; sex=0)
        plan = Plan(25, 50, 200)
        fixed = ["grt"]
        dF = 0.011
        F0 = 0.027

        foo, bar, tstDir = "rand", "dsblup", "rst/onetrt"
        xy = "$tstDir/$bar.xy"
        cp("$tstDir/$foo.xy", xy, force=true)
        ped = deserialize("$tstDir/$foo.ped")
        G = zeros(nrow(ped), nrow(ped))
        read!("$tstDir/$foo.irm", G)

        ids = view(ped, ped.grt .== ped.grt[end], :id) # ID of the last generation
        phenotype!(ids, ped, milk)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, milk)
        g22 = G[ids, ids]
        mid = nrow(ped)
        ng = Select(ids, ped, g22, milk, plan.noff, dF, 1; F0=F0, ocs=TM2024)
        #reproduce!(ng, ped, xy, lmp, milk)
        #G = xirm(G, xy, lmp.chip, mid, nrow(ped))
    elseif op == 3
    else
    end
end
