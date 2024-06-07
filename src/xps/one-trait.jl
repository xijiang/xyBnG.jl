function one_trait()
    # Scenarios
    baseDir, tstDir = "rst/base", "rst/onetrt"
    cattle = Cattle("cattle", 5_000)
    milk = Trait("milk", 0.25, 10000; sex = 0)
    nid, nchp, nref = 200, 50_000, 10_000
    nrng, nsel = 5, 10
    plan = Plan(25, 50, 200) # default :hierarchical
    fixed = ["grt"]

    title = md"""# Simulation: Selection on one trait"""
    tprintln(title)
    synopsis = md"""## Synpsis:
    In this play, a founder cattle population is sampled from a base population.
    It is then randomly selected on milk for 5 generations.
    The DOSc methods is then tested with the random selection results.
    """
    tprintln(synopsis, "\n")

    @info "Prepare a founder population"
    println("\n")

    isfile("$baseDir/desc.txt") || sim_base(cattle, baseDir)
    species = readline("$baseDir/desc.txt")
    if !isfile("$tstDir/founder.xy")
        sample_founder(baseDir, tstDir, nid, nchp, nref, milk)
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
    tprintln(act)
    ped = deserialize("$tstDir/founder.ped")
    foo, bar = "founder", "rand"
    cp("$tstDir/$foo.xy", "$tstDir/$bar.xy", force=true)
    xy = "$tstDir/$bar.xy"
    
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

    act = md"""
    ## Act II: Directional selection on milk
    ### Scene 1: AABLUP method
    """
    tprintln(act)
    foo, bar = "rand", "aablup"
    cp("$tstDir/$foo.xy", "$tstDir/$bar.xy", force=true)
    xy = "$tstDir/$bar.xy"
    ped = deserialize("$tstDir/rand.ped")
    @info "  - Operating on generation:"
    for igrt in 1:nsel
        print(" $igrt")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, milk)
        G = nrm(ped)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, milk)
        ng = Select(ids, plan, ped, milk)
        reproduce!(ng, ped, xy, lmp, milk)
    end
    println("\n")
    serialize("$tstDir/$bar.ped", ped)

    act = md"""
    ### Scene 2: DOSc method
    """
    tprintln(act)
    #@goto debug_skip
    #@label debug_skip
end
