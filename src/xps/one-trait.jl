function one_trait(;
                    baseDir = "rst/base",
                    tstDir = "rst/onetrt",
                    cattle = Cattle("cattle", 5_000),
                    milk = Trait("milk", 0.25, 10000; sex = 0),
                    nid = 200,
                    nchp = 50_000,
                    nref = 10_000,
                    nrng = 5,
                    nsel = 20,
                    plan = Plan(25, 50, 200),
                    fixed = ["grt"],
                    dF = 0.011,)
    title = md"""# Simulation: Selection on one trait"""
    tprintln(title, "\n")
    synopsis = md"""## Synpsis:
    In this play, a founder cattle population is sampled from a base population.
    It is then randomly selected on milk for 5 generations.
    The DOSc methods is then tested with the random selection results.
    """
    tprintln(synopsis, "\n")

    @info "Prepare a founder population"
    println("\n")

    isfile("$baseDir/desc.txt") || ts_base(cattle, baseDir)
    species = readline("$baseDir/desc.txt")
    if !isfile("$tstDir/founder.xy")
        sample_ts(baseDir, tstDir, nid, nchp, nref, milk)
        for ext in ["ped", "lmp"]
            mv("$tstDir/$species.$ext", "$tstDir/founder.$ext", force=true)
        end
        uniq("$tstDir/$species.xy", "$tstDir/founder.xy")
        rm("$tstDir/$species.xy", force=true)
    end
    lmp = deserialize("$tstDir/founder.lmp")

    act = md"""
    ## Act I: Random selection on milk
    """
    tprintln(act, "\n")
    ped = deserialize("$tstDir/founder.ped")
    foo, bar = "founder", "rand"
    xy = "$tstDir/$bar.xy"
    cp("$tstDir/$foo.xy", xy, force=true)
    
    @info "  - Operating on generation:"
    for igrt in 1:nrng
        print(" $igrt")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, milk)
        Predict!(ids, ped, milk)
        ng = Select(ids, plan, ped, milk)
        reproduce!(ng, ped, xy, lmp, milk)
    end
    println("\n")
    serialize("$tstDir/$bar.ped", ped)
    @info "  - Calculating IBD relationship matrix"
    println("\n")
    G = irm(xy, lmp.chip, ped.id)
    write("$tstDir/$bar.irm", G)
    F0 = 0.027

    act = md"""
    ## Act II: Directional selection on milk
    
    ### Scene 1: AABLUP method
    """
    tprintln(act, "\n")
    foo, bar = "rand", "aablup"
    xy = "$tstDir/$bar.xy"
    cp("$tstDir/$foo.xy", xy, force=true)
    ped = deserialize("$tstDir/$foo.ped")
    @info "  - Operating on generation:"
    for igrt in 1:nsel
        print(" $igrt")
        ids = view(ped, ped.grt .== ped.grt[end], :id) # ID of the last generation
        phenotype!(ids, ped, milk)
        G = nrm(ped)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, milk)
        g22 = G[ids, ids]
        ng = Select(ids, ped, g22, milk, plan.noff, dF, igrt; F0 = F0)
        reproduce!(ng, ped, xy, lmp, milk)
    end
    println("\n")
    serialize("$tstDir/$bar.ped", ped)

    act = md"""
    ### Scene 2: IIBLUP method
    """
    tprintln(act, "\n")
    foo, bar = "rand", "iiblup"
    xy = "$tstDir/$bar.xy"
    cp("$tstDir/$foo.xy", xy, force=true)
    ped = deserialize("$tstDir/$foo.ped")
    G = zeros(nrow(ped), nrow(ped))
    read!("$tstDir/$foo.irm", G)
    @info "  - Operating on generation:"
    for igrt in 1:nsel
        print(" $igrt")
        ids = view(ped, ped.grt .== ped.grt[end], :id) # ID of the last generation
        phenotype!(ids, ped, milk)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, milk)
        g22 = G[ids, ids]
        mid = nrow(ped)
        ng = Select(ids, ped, g22, milk, plan.noff, dF, igrt; F0 = F0)
        reproduce!(ng, ped, xy, lmp, milk)
        # update the IBD relationship matrix
        G = xirm(G, xy, lmp.chip, mid, nrow(ped))
    end
    println("\n")
    serialize("$tstDir/$bar.ped", ped)

    #=
    act = md"""
    ### Scene 3: DOSc method
    """
    tprintln(act, "\n")
    foo, bar = "rand", "dsblup"
    xy = "$tstDir/$bar.xy"
    cp("$tstDir/$foo.xy", xy, force=true)
    ped = deserialize("$tstDir/$foo.ped")
    G = zeros(nrow(ped), nrow(ped))
    read!("$tstDir/$foo.irm", G)
    @info "  - Operating on generation:"
    for igrt in 1:nsel
        print(" $igrt")
        ids = view(ped, ped.grt .== ped.grt[end], :id) # ID of the last generation
        phenotype!(ids, ped, milk)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, milk)
        g22 = G[ids, ids]
        mid = nrow(ped)
        @debug mid
        ng = Select(ids, ped, g22, milk, plan.noff, dF, igrt;
                    F0 = F0, ocs = TM2024)
        reproduce!(ng, ped, xy, lmp, milk)
        G = xirm(G, xy, lmp.chip, mid, nrow(ped))
    end
    =#
end
