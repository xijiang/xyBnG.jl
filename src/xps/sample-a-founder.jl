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
    for ext in ["ped", "lmp", "xy"]
        cp("$fdrDir/cattle.$ext", "$fdrDir/founder.$ext")
    end
end
