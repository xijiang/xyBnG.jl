# This are generally some example selection schemes
# They should be easy to setup using the graph expression of breeding programs
# !!! All these programs are selection the last generation of the pedigree
# the first letter of the function is the matrix used for constraint
# the letter next to `ocs` is the matrix used for EBV
# Parameters of these functions are:
# - `test` is the directory where the simulation is stored.
# - `foo` is the prefix of the pedigree and genotype files to be selected.
# - `bar` is the prefix of the selected pedigree and genotype files to be saved.
# - `lmp` is the linkage map.
# - `ngn` is the number of generations to be selected.
# - `trait` is the trait to be selected.
# - `fixed` are the fixed effects.
# - `noff` is the number of offspring.
# - `dF` is the ΔF of each generation.
# - `F0` is the inbreeding coefficient of the initial population.
# The selected parents are randomly selected on contibution weights to mate.

"""
    savepar(dic, par)
Save the parameters of the current simulation into file `par`.
"""
function savepar(dic, par)
    for (key, val) in dic
        @show key, val
    end
end

"""
    randbrd(test, foo, bar, lmp, ngn, trait, plan; ibd = false)
Random selection on `foo`.xy, `foo`.ped in directory `test` for `ngn`
generations. This is to select a single trait `trait` using a selection plan
`plan`. SNP linkage information are in DataFrame `lmp`. The results are saved in
`bar`.xy, `bar`.ped in directory `test`. If `ibd` is `true`, the IBD
relationship matrix is calculated and saved in `bar`.irm.
"""
function randbrd(test, foo, bar, lmp, ngn, trait, plan; ibd = false)
    @info "  - Random selection for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force=true)
    for ign in 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        Predict!(ids, ped, trait)
        ng = Select(ids, plan, ped, trait)
        reproduce!(ng, ped, xy, lmp, trait)
    end
    println()
    serialize("$test/$bar.ped", ped)
    if ibd
        @info "  - Calculating IBD relationship matrix"
        G = irm(xy, lmp.chip, 1:nrow(ped))
        write("$test/$bar.irm", G)
    end
end

"""
    aaocs(test, foo, bar, lmp, ngn, trait, fixed, noff, dF, F0)
Optimal contribution selection with `A` for both EBV and constraint on `foo`.xy
and `foo`.ped in directory `test` for `ngn` generations. SNP linkage information
are in DataFrame `lmp`. The results are saved in `bar`.xy, `bar`.ped in
directory `test`. The selection is on a single trait `trait` with fixed effects
`fixed`, which is a column name vector in pedigree DataFrame. The `noff` pair
of parents are sampled according to their contribution weights. The constraint
ΔF is `dF`. `F0` is the inbreeding coefficient of the `foo` population.


See also [`randbrd`](@ref), [`iiocs`](@ref), [`iidos`](@ref), [`ggocs`](@ref),
[`agocs`](@ref), [`igocs`](@ref).
"""
function aaocs(test, foo, bar, lmp, ngn, trait, fixed, noff, dF, F0)
    @info "  - Directional selection AABLUP for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force=true)
    for ign in 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        G = nrm(ped)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        g22 = G[ids, ids]
        ng = Select(ids, ped, g22, trait, noff, dF, ign; F0=F0)
        reproduce!(ng, ped, xy, lmp, trait)
    end
    println()
    serialize("$test/$bar.ped", ped)
end

"""
    iiocs(test, foo, bar, lmp, ngn, trait, fixed, noff, dF, F0)
Optimal contribution selection with `IBD` relationship matrix for both EBV and
constraint on `foo`.xy and `foo`.ped in directory `test` for `ngn` generations.
SNP linkage information are in DataFrame `lmp`. The results are saved in
`bar`.xy, `bar`.ped in directory `test`. The selection is on a single trait
`trait` with fixed effects `fixed`, which is a column name vector in pedigree
DataFrame. The `noff` pair of parents are sampled according to their
contribution weights. The constraint ΔF is `dF`. `F0` is the inbreeding
coefficient of the `foo` population.

This function uses the TM1997 algorithm for OCS.

See also [`randbrd`](@ref), [`aaocs`](@ref), [`iidos`](@ref), [`ggocs`](@ref),
[`agocs`](@ref), [`igocs`](@ref).
"""
function iiocs(test, foo, bar, lmp, ngn, trait, fixed, noff, dF, F0)
    @info "  - Directional selection IIBLUP for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force=true)
    G = zeros(nrow(ped), nrow(ped))
    read!("$test/$foo.irm", G)
    for ign in 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        g22 = G[ids, ids]
        mid = nrow(ped)
        ng = Select(ids, ped, g22, trait, noff, dF, ign; F0=F0)
        reproduce!(ng, ped, xy, lmp, trait)
        G = xirm(G, xy, lmp.chip, mid, nrow(ped))
    end
    println()
    serialize("$test/$bar.ped", ped)
end

"""
    iidos(test, foo, bar, lmp, ngn, trait, fixed, noff, dF, F0)
Optimal contribution selection with `IBD` relationship matrix for both EBV and
constraint on `foo`.xy and `foo`.ped in directory `test` for `ngn` generations.
SNP linkage information are in DataFrame `lmp`. The results are saved in
`bar`.xy, `bar`.ped in directory `test`. The selection is on a single trait
`trait` with fixed effects `fixed`, which is a column name vector in pedigree
DataFrame. The `noff` pair of parents are sampled according to their
contribution weights. The constraint ΔF is `dF`. `F0` is the inbreeding
coefficient of the `foo` population.

This function uses the TM2024 algorithm for OCS.

See also [`randbrd`](@ref), [`aaocs`](@ref), [`iiocs`](@ref), [`ggocs`](@ref),
[`agocs`](@ref), [`igocs`](@ref).
"""
function iidos(test, foo, bar, lmp, ngn, trait, fixed, noff, dF, F0)
    @info "  - Directional selection DOSc for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force=true)
    G = zeros(nrow(ped), nrow(ped))
    read!("$test/$foo.irm", G)
    for ign in 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        g22 = G[ids, ids]
        mid = nrow(ped)
        ng = Select(ids, ped, g22, trait, noff, dF, ign; F0=F0, ocs=TM2024)
        reproduce!(ng, ped, xy, lmp, trait)
        G = xirm(G, xy, lmp.chip, mid, nrow(ped))
    end
    println()
    serialize("$test/$bar.ped", ped)
end

"""
    ggocs(test, foo, bar, lmp, ngn, trait, fixed, noff, dF, F0)
Optimal contribution selection with `G` relationship matrix for both EBV and
constraint on `foo`.xy and `foo`.ped in directory `test` for `ngn` generations.
SNP linkage information are in DataFrame `lmp`. The results are saved in
`bar`.xy, `bar`.ped in directory `test`. The selection is on a single trait
`trait` with fixed effects `fixed`, which is a column name vector in pedigree
DataFrame. The `noff` pair of parents are sampled according to their
contribution weights. The constraint ΔF is `dF`. `F0` is the inbreeding
coefficient of the `foo` population.

This function uses the TM1997 algorithm for OCS.

See also [`randbrd`](@ref), [`aaocs`](@ref), [`iiocs`](@ref), [`iiocs`](@ref),
[`agocs`](@ref), [`igocs`](@ref).
"""
function ggocs(test, foo, bar, lmp, ngn, trait, fixed, noff, dF, F0)
    @info "  - Directional selection GGBLUP for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force=true)
    for ign in 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        G = grm(xy, lmp.chip)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        g22 = G[ids, ids]
        ng = Select(ids, ped, g22, trait, noff, dF, ign; F0=F0, ong=true)
        reproduce!(ng, ped, xy, lmp, trait)
    end
    println()
    serialize("$test/$bar.ped", ped)
end

"""
    agocs(test, foo, bar, lmp, ngn, trait, fixed, noff, dF, F0)
Optimal contribution selection with `G` relationship matrix for EBV, `A`
relationship matrix for constraint on `foo`.xy and `foo`.ped in directory `test`
for `ngn` generations. SNP linkage information are in DataFrame `lmp`. The
results are saved in `bar`.xy, `bar`.ped in directory `test`. The selection is
on a single trait `trait` with fixed effects `fixed`, which is a column name
vector in pedigree DataFrame. The `noff` pair of parents are sampled according
to their contribution weights. The constraint ΔF is `dF`. `F0` is the inbreeding
coefficient of the `foo` population.

This function uses the TM1997 algorithm for OCS.

See also [`randbrd`](@ref), [`aaocs`](@ref), [`iiocs`](@ref), [`iiocs`](@ref),
[`ggocs`](@ref), [`igocs`](@ref).
"""
function agocs(test, foo, bar, lmp, ngn, trait, fixed, noff, dF, F0)
    @info "  - Directional selection AGBLUP for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force=true)
    for ign in 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        G = grm(xy, lmp.chip)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        G = nrm(ped)
        g22 = G[ids, ids]
        ng = Select(ids, ped, g22, trait, noff, dF, ign; F0=F0)
        reproduce!(ng, ped, xy, lmp, trait)
    end
    println()
    serialize("$test/$bar.ped", ped)
end

"""
    igocs(test, foo, bar, lmp, ngn, trait, fixed, noff, dF, F0)
Optimal contribution selection with `G` relationship matrix for EBV, `IBD`
relationship matrix for constraint on `foo`.xy and `foo`.ped in directory `test`
for `ngn` generations. SNP linkage information are in DataFrame `lmp`. The
results are saved in `bar`.xy, `bar`.ped in directory `test`. The selection is
on a single trait `trait` with fixed effects `fixed`, which is a column name
vector in pedigree DataFrame. The `noff` pair of parents are sampled according
to their contribution weights. The constraint ΔF is `dF`. `F0` is the inbreeding
coefficient of the `foo` population.

This function uses the TM1997 algorithm for OCS.

See also [`randbrd`](@ref), [`aaocs`](@ref), [`iiocs`](@ref), [`iiocs`](@ref),
[`ggocs`](@ref), [`agocs`](@ref).
"""
function igocs(test, foo, bar, lmp, ngn, trait, fixed, noff, dF, F0)
    @info "  - Directional selection IGBLUP for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force=true)
    for ign in 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        G = grm(xy, lmp.chip)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        mid = nrow(ped)
        g22 = irm(xy, lmp.chip, mid+1-length(ids):mid)
        ng = Select(ids, ped, g22, trait, noff, dF, ign; F0=F0)
        reproduce!(ng, ped, xy, lmp, trait)
    end
    println()
    serialize("$test/$bar.ped", ped)
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
    cp("$test/$foo.xy", xy, force=true)
    for ign in 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        G = grm(xy, lmp.chip)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        ng = Select(ids, plan, ped, trait)
        reproduce!(ng, ped, xy, lmp, trait)
    end
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
    cp("$test/$foo.xy", xy, force=true)
    for ign in 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        G = nrm(ped)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        ng = Select(ids, plan, ped, trait)
        reproduce!(ng, ped, xy, lmp, trait)
    end
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
    cp("$test/$foo.xy", xy, force=true)
    G = zeros(nrow(ped), nrow(ped))
    read!("$test/$foo.irm", G)
    for ign in 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        g22 = G[ids, ids]
        mid = nrow(ped)
        ng = Select(ids, plan, ped, trait)
        reproduce!(ng, ped, xy, lmp, trait)
        G = xirm(G, xy, lmp.chip, mid, nrow(ped))
    end
    println()
    serialize("$test/$bar.ped", ped)
end
