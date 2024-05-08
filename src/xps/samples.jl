"""
    function sim_cattle_base()

This sample function is to simulate a cattle base population, using `stdpopsim`.
"""
function sim_cattle_base()
    cattle = Cattle("cattle", 5_000)
    sim_base(cattle, "rst/base")
end

"""
    sample_a_founder()

This sample function is to sample a founder population of 100 ID, 10,000 chip
SNPs, 5,000 reference SNPs and two traits, milk and growth, from result of
`base()` into a directory "sandbox".
"""
function sample_a_founder()
    milk   = Trait("milk",   0.3, 10000; sex = 0, vd = 0.1)
    growth = Trait("growth", 0.5, 10000; sex = 2, vd = 0.2)
    nid, nchp, nref = 200, 50_000, 10_000
    baseDir, fdrDir = "rst/base", "rst/cattle"
    sample_founder(baseDir, fdrDir, nid, nchp, nref, milk, growth)
end

#=
"""
    rand_sel_one(bar)

Random selection on one trait.  The trait's EBV are randomly assigned.
"""
function rand_sel(bar)
    milk   = Trait("milk",   0.3, 1000, vd = 0.1)
    growth = Trait("growth", 0.5, 1000, vd = 0.2)
    dir = "rst/sandbox"
    lmp = deserialize("$dir/$bar-map.ser")
    npa, nma, noff, ngrt = 10, 20, 100, 5

    ped = deserialize("$dir/$bar-ped.ser")
    fxy = "$dir/$bar-milk.xy"
    cp("$dir/$bar-fdr.xy", fxy, force = true)
    for _ in 1:ngrt
        predict!(ped, milk)
        Select!(ped, (npa, nma, noff), milk; mate = :balanced)
        reproduce!(ped, fxy, lmp, milk, growth)
    end
    serialize("$dir/$bar-milk.ped", ped)
    ped
end

"""
    rand_sel_two(bar)

Random selection on both traits. create new files for selection. One can of
course do selection from the milk selection results above.
"""
function rand_sel_two(bar)
    milk   = Trait("milk",   0.3, 1000, vd = 0.1)
    growth = Trait("growth", 0.5, 1000, vd = 0.2)
    dir = "rst/sandbox"
    lmp = deserialize("$dir/$bar-map.ser")
    npa, nma, noff, ngrt = 10, 20, 100, 5

    ped = deserialize("$dir/$bar-ped.ser")
    fxy = "$dir/$bar-both.xy"
    cp("$dir/$bar-fdr.xy", fxy, force = true)
    wgt = Dict(:milk => 0.5, :growth => 0.5)
    for _ in 1:ngrt
        predict!(ped, milk, growth)
        Select!(ped, (npa, nma, noff), wgt; mate = :balanced)
        reproduce!(ped, fxy, lmp, milk, growth)
    end
    serialize("$dir/$bar-both.ped", ped)
    ped
end
=#