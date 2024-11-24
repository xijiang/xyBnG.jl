"""
    testScheme()
This function is to run a very small example of the package `xyBnG`. It is
for my laptop to run in a reasonable time.
"""
function tstScheme(test, scheme)
    foo, bar = "mini-rand", string(scheme)
    @info "Test scheme: $bar"
    lmp = deserialize("$test/mini-founder.lmp")
    trait = Trait("growth", 0.25, 1_000)
    nsel, fixed = 5, ["grt"]
    plan = Plan(5, 10, 40)
    dF, F0, ε = 0.03, 0.14177145996093743, 1e-6
    if scheme ∈ (iiocs, ggocs, igocs, hgocs, hhocs)
        scheme(test, foo, bar, lmp, nsel, trait, fixed, plan, dF, F0; ε = ε)
    elseif scheme ∈ [aaocs]
        scheme(test, foo, bar, lmp, nsel, trait, fixed, plan, dF, F0)
    elseif scheme ∈ [gblup, iblup]
        scheme(test, foo, bar, lmp, nsel, trait, fixed, plan; ε = ε)
    elseif scheme ∈ [ablup]
        scheme(test, foo, bar, lmp, nsel, trait, fixed, plan)
    else
        error("Scheme $scheme is not defined.")
    end
end

"""
    minifdr(base, tdir)
Sample a mini founder population (5 x 10 → 40, each generation) from `base` with
3k chip SNP, 1k reference SNP and 1k QTL for a growth trait of h2 =0.25. It then
run 5 generations of random selection. The founder population is saved in
`tdir`.

## Example
```julia
xyBnG.xps.minifdr("/home/xijiang/workspace/base/tskit", "/home/xijiang/workspace/scheme/")
````
"""
function minifdr(base, tdir)
    species = Cattle(5_000)
    chkbase(base, species) # Prepare/verify a base population
    fxy, fmp, maf = "$base/$(species.name).xy", "$base/$(species.name).lmp", 0.0
    plan = Plan(5, 10, 40)
    trait = Trait("growth", 0.25, 1_000)
    nchp, nref, nrng = 3_000, 1_000, 5
    initPop(fxy, fmp, tdir, plan, maf, nchp, nref, nrng, trait, "mini")
end
