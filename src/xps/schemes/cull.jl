"""
    randbrd(test, foo, bar, lmp, ngn, trait, plan)
Random selection on `foo`.xy, `foo`.ped in directory `test` for `ngn`
generations. This is to select a single trait `trait` using a selection plan
`plan`. SNP linkage information are in DataFrame `lmp`. The results are saved in
`bar`.xy, `bar`.ped in directory `test`. If `ibd` is `true`, the IBD
relationship matrix is calculated and saved in `bar`.irm.
"""
function randbrd(test, foo, bar, lmp, ngn, trait, plan)
    @info "  - Random selection for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force = true)
    for ign = 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        Predict!(ids, ped, trait)
        ng = Select(ids, plan, ped, trait)
        reproduce!(ng, ped, xy, lmp, trait)
    end
    println()
    serialize("$test/$bar.ped", ped)
    δF = 1 / 8plan.npa + 1 / 8plan.nma  # Falconer 1996, pp. 67
    n = ngn > 2 ? ngn - 1 : 0
    1 - (1 - δF)^n # F0: expected inbreeding after random selection
end

"""
    gblup(test, foo, bar, lmp, ngn, trait, fixed, plan)
Directional selection with `G` relationship matrix for EBV on `foo`.xy and
`foo`.ped in directory `test` for `ngn` generations. SNP linkage information are
in DataFrame `lmp`. The results are saved in `bar`.xy, `bar`.ped in directory
`test`. The selection is on a single trait `trait` with fixed effects `fixed`,
which is a column name vector in pedigree DataFrame. The selection is based on
the selection plan `plan`.

See also [`ablup`](@ref), [`iblup`](@ref).
"""
function gblup(test, foo, bar, lmp, ngn, trait, fixed, plan)
    @info "  - Directional selection GBLUP for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force = true)
    for ign = 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        G = grm(xy, lmp.chip, lmp.frq)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        ng = Select(ids, plan, ped, trait)
        reproduce!(ng, ped, xy, lmp, trait)
    end
    println()
    serialize("$test/$bar.ped", ped)
end

"""
    ablup(test, foo, bar, lmp, ngn, trait, fixed, plan)
Directional selection with `A` relationship matrix for EBV on `foo`.xy and
`foo`.ped in directory `test` for `ngn` generations. SNP linkage information are
in DataFrame `lmp`. The results are saved in `bar`.xy, `bar`.ped in directory
`test`. The selection is on a single trait `trait` with fixed effects `fixed`,
which is a column name vector in pedigree DataFrame. The selection is based on
the selection plan `plan`.

See also [`gblup`](@ref), [`iblup`](@ref).
"""
function ablup(test, foo, bar, lmp, ngn, trait, fixed, plan)
    @info "  - Directional selection ABLUP for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force = true)
    for ign = 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        G = nrm(ped)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        ng = Select(ids, plan, ped, trait)
        reproduce!(ng, ped, xy, lmp, trait)
    end
    println()
    serialize("$test/$bar.ped", ped)
end

"""
    iblup(test, foo, bar, lmp, ngn, trait, fixed, plan)
Directional selection with `IBD` relationship matrix for EBV on `foo`.xy and
`foo`.ped in directory `test` for `ngn` generations. SNP linkage information are
in DataFrame `lmp`. The results are saved in `bar`.xy, `bar`.ped in directory
`test`. The selection is on a single trait `trait` with fixed effects `fixed`,
which is a column name vector in pedigree DataFrame. The selection is based on
the selection plan `plan`.

See also [`gblup`](@ref), [`ablup`](@ref).
"""
function iblup(test, foo, bar, lmp, ngn, trait, fixed, plan)
    @info "  - Directional selection IBLUP for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force = true)
    G = nothing
    if isfile("$test/$foo.irm")
        G = zeros(size(ped, 1), size(ped, 1))
        read!("$test/$foo.irm", G)
    else
        @info "  - Calculating IBD relationship matrix"
        G = irm(xy, lmp.chip, 1:size(ped, 1))
        write("$test/$foo.irm", G)
    end
    for ign = 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        mid = size(ped, 1)
        ng = Select(ids, plan, ped, trait)
        reproduce!(ng, ped, xy, lmp, trait)
        G = xirm(G, xy, lmp.chip, mid, size(ped, 1))
    end
    println()
    serialize("$test/$bar.ped", ped)
end
