"""
    initPop(
        fxy::AbstractString,
        fmp::AbstractString,
        test::AbstractString,
        plan::Plan,
        maf::Float64,
        nchp::Int,
        nref::Int,
        nrng::Int,
        trait::Trait,
        tag::AbstractString,
        cibd::Bool,
    )
Sample a founder population and run a breeding scheme with random selection.
"""
function initPop(
    fxy::AbstractString,
    fmp::AbstractString,
    test::AbstractString,
    plan::Plan,
    maf::Float64,
    nchp::Int,
    nref::Int,
    nrng::Int,
    trait::Trait,
    tag::AbstractString,
    cibd::Bool,
)
    sample_xy(fxy, fmp, test, plan.noff, maf, nchp, nref, trait)
    mv("$test/founder.xy", "$test/$tag-snp.xy", force = true)
    mv("$test/founder.lmp", "$test/$tag-founder.lmp", force = true)
    mv("$test/founder.ped", "$test/$tag-founder.ped", force = true)
    uniq("$test/$tag-snp.xy", "$test/$tag-founder.xy")
    lmp = deserialize("$test/$tag-founder.lmp")
    randbrd(test, "$tag-founder", "$tag-rand", lmp, nrng, trait, plan; ibd = cibd)
    lmp
end
