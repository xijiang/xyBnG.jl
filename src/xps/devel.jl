using xyBnG.Breeding
using Statistics
import StatsBase: sample, Weights

function three(;
    test = "sim/debug",
    trait = Trait("milk", 0.25, 10000; sex = 0),
    nrng = 5,
    nsel = 10,
    plan = Plan(25, 50, 200),
    fixed = ["grt"],
    dF = 0.011,
    F0 = 0.027,
)
    for irpt in 1:10
        # Sample a founder population
        tskit = "sim/tskit"
        fxy, fmp, maf = "$tskit/BosTau.xy", "$tskit/BosTau.lmp", 0.0
        sample_xy(fxy, fmp, test, 200, maf, 50000, 10000, trait)
        mv("$test/founder.xy", "$test/snp.xy", force = true)
        uniq("$test/snp.xy", "$test/founder.xy")
        
        ## Test the random selection
        lmp = deserialize("$test/founder.lmp")
        randbrd(test, "founder", "rand", lmp, nrng, trait, plan)

        #schemes = (aaocs, iiocs, iidos, ggocs, agocs, igocs, gblup, ablup, iblup)
        foo, bar = "rand", "ggocs"
        ggocs(test, foo, bar, lmp, nsel, trait, fixed, plan.noff, dF, F0)
        ped = deserialize("$test/$bar.ped")
        df = filter(row -> 5 < row.grt <15, ped)
        @info "Repeat: $irpt"
        for x in groupby(df, :grt)
            println(x.grt[end], ' ', 
                    length(unique(x.sire)), '\t', 
                    length(unique(x.dam)), '\t',)
                    #round(var(x.ebv_milk[x.sex .== 1]), digits=3), '\t',
                    #round(var(x.ebv_milk[x.sex .== 0]), digits=3))
        end
    end
end
