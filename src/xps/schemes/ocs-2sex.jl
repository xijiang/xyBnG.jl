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
    aaocs(test, foo, bar, lmp, ngn, trait, fixed, plan, dF, F0)
Optimal contribution selection with `A` for both EBV and constraint on `foo`.xy
and `foo`.ped in directory `test` for `ngn` generations. SNP linkage information
are in DataFrame `lmp`. The results are saved in `bar`.xy, `bar`.ped in
directory `test`. The selection is on a single trait `trait` with fixed effects
`fixed`, which is a column name vector in pedigree DataFrame. Parents are
sampled according to `plan`. The constraint ΔF is `dF`. `F0` is the inbreeding
coefficient of the `foo` population.


See also [`randbrd`](@ref), [`iiocs`](@ref), [`iidos`](@ref), [`ggocs`](@ref),
[`agocs`](@ref), [`igocs`](@ref).
"""
function aaocs(test, foo, bar, lmp, ngn, trait, fixed, plan, dF, F0)
    @info "  - Directional selection AABLUP for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force = true)
    for ign = 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        G = nrm(ped)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        g22 = G[ids, ids]
        ng = Select(ids, plan, ped, g22, trait, dF, ign; F0 = F0)
        reproduce!(ng, ped, xy, lmp, trait)
    end
    println()
    serialize("$test/$bar.ped", ped)
end

"""
    iiocs(test, foo, bar, lmp, ngn, trait, fixed, plan, dF, F0; ε = 1e-6)
Optimal contribution selection with `IBD` relationship matrix for both EBV and
constraint on `foo`.xy and `foo`.ped in directory `test` for `ngn` generations.
SNP linkage information are in DataFrame `lmp`. The results are saved in
`bar`.xy, `bar`.ped in directory `test`. The selection is on a single trait
`trait` with fixed effects `fixed`, which is a column name vector in pedigree
DataFrame. Parents are sampled according to `plan`. The constraint ΔF is `dF`.
`F0` is the inbreeding coefficient of the `foo` population.

This function uses the TM1997 algorithm for OCS.

See also [`randbrd`](@ref), [`aaocs`](@ref), [`iidos`](@ref), [`ggocs`](@ref),
[`agocs`](@ref), [`igocs`](@ref).
"""
function iiocs(test, foo, bar, lmp, ngn, trait, fixed, plan, dF, F0; ε = 1e-6)
    @info "  - Directional selection IIOCS for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force = true)
    G = fileIRM("$test/$foo.irm", xy, lmp.chip, 1:size(ped, 1); ε = ε)
    for ign = 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        g22 = G[ids, ids]
        mid = size(ped, 1)
        ng = Select(ids, plan, ped, g22, trait, dF, ign; F0 = F0)
        reproduce!(ng, ped, xy, lmp, trait)
        G = xirm(G, xy, lmp.chip, mid, size(ped, 1))
        for i = mid+1:size(G, 1)
            G[i, i] += ε
        end
    end
    println()
    serialize("$test/$bar.ped", ped)
end

"""
    riocs(test, foo, bar, lmp, ngn, trait, fixed, plan, dF, F0; ε = 1e-6)
Optimal contribution selection with `IBD` relationship matrix for both EBV and
constraint on `foo`.xy and `foo`.ped in directory `test` for `ngn` generations.
For constraint, the IRM is from the reference SNP loci. The EBV are estimated
with IRM from chip SNPs. SNP linkage information are in DataFrame `lmp`. The
results are saved in `bar`.xy, `bar`.ped in directory `test`. The selection is
on a single trait `trait` with fixed effects `fixed`, which is a column name
vector in pedigree DataFrame. Parents are sampled according to `plan`. The
constraint ΔF is `dF`. `F0` is the inbreeding coefficient of the `foo`
population.

This function uses the TM1997 algorithm for OCS.

See also [`randbrd`](@ref), [`aaocs`](@ref), [`iidos`](@ref), [`ggocs`](@ref),
[`agocs`](@ref), [`igocs`](@ref).
"""
function riocs(test, foo, bar, lmp, ngn, trait, fixed, plan, dF, F0; ε = 1e-6)
    @info "  - Directional selection IIOCS for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force = true)
    G = fileIRM("$test/$foo.irm", xy, lmp.dark, 1:size(ped, 1); ε = ε)
    for ign = 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        g22 = G[ids, ids]
        mid = size(ped, 1)
        ng = Select(ids, plan, ped, g22, trait, dF, ign; F0 = F0)
        reproduce!(ng, ped, xy, lmp, trait)
        G = xirm(G, xy, lmp.dark, mid, size(ped, 1))
        for i = mid+1:size(G, 1)
            G[i, i] += ε
        end
    end
    println()
    serialize("$test/$bar.ped", ped)
end

"""
    iidos(test, foo, bar, lmp, ngn, trait, fixed, plan, dF, F0)
Optimal contribution selection with `IBD` relationship matrix for both EBV and
constraint on `foo`.xy and `foo`.ped in directory `test` for `ngn` generations.
SNP linkage information are in DataFrame `lmp`. The results are saved in
`bar`.xy, `bar`.ped in directory `test`. The selection is on a single trait
`trait` with fixed effects `fixed`, which is a column name vector in pedigree
DataFrame. Parents are sampled according to `plan`. The constraint ΔF is `dF`.
`F0` is the inbreeding coefficient of the `foo` population.

This function uses the TM2024 algorithm for OCS.

See also [`randbrd`](@ref), [`aaocs`](@ref), [`iiocs`](@ref), [`ggocs`](@ref),
[`agocs`](@ref), [`igocs`](@ref).
"""
function iidos(test, foo, bar, lmp, ngn, trait, fixed, plan, dF, F0)
    @info "  - Directional selection DOSc for $ngn generations"
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
        g22 = G[ids, ids]
        mid = size(ped, 1)
        ng = Select(ids, plan, ped, g22, trait, dF, ign; F0 = F0, ocs = TM2024)
        reproduce!(ng, ped, xy, lmp, trait)
        G = xirm(G, xy, lmp.chip, mid, size(ped, 1))
    end
    println()
    serialize("$test/$bar.ped", ped)
end

#=
"""
    tgocs(test, foo, bar, lmp, ngn, trait, fixed, noff, dF, F0)
Optimal contribution selection with `G` relationship matrix for both EBV and
constraint on `foo`.xy and `foo`.ped in directory `test` for `ngn` generations.
SNP linkage information are in DataFrame `lmp`. The results are saved in
`bar`.xy, `bar`.ped in directory `test`. The selection is on a single trait
`trait` with fixed effects `fixed`, which is a column name vector in pedigree
DataFrame. Parent are sampled according to `plan`. The constraint ΔF is `dF`.
`F0` is the inbreeding coefficient of the `foo` population.

This function uses the TM1997 algorithm for OCS. As the constraint matrix is
using GRM, which has already considered the allele frequency changes, hence
option `ong` is set to `true`.

See also [`randbrd`](@ref), [`aaocs`](@ref), [`iiocs`](@ref), [`iiocs`](@ref),
[`agocs`](@ref), [`igocs`](@ref).
"""
function tgocs(test, foo, bar, lmp, ngn, trait, fixed, plan, dF, F0)
    @info "  - Directional selection TGBLUP for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force = true)
    for ign = 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        G = grm(xy, lmp.chip, lmp.frq)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        g22 = G[ids, ids]
        ng = Select(ids, plan, ped, g22, trait, dF, ign; F0 = F0, ong = true)
        reproduce!(ng, ped, xy, lmp, trait)
    end
    println()
    serialize("$test/$bar.ped", ped)
end
=#

"""
    ggocs(test, foo, bar, lmp, ngn, trait, fixed, plan, dF, F0; ε = 1e-6)
Optimal contribution selection with `G` relationship matrix for both EBV and
constraint on `foo`.xy and `foo`.ped in directory `test` for `ngn` generations.
SNP linkage information are in DataFrame `lmp`. The results are saved in
`bar`.xy, `bar`.ped in directory `test`. The selection is on a single trait
`trait` with fixed effects `fixed`, which is a column name vector in pedigree
DataFrame. Parents are sampled according to `plan`. The constraint ΔF is `dF`.
`F0` is the inbreeding coefficient of the `foo` population.

This function uses the TM1997 algorithm for OCS. As the constraint matrix is
using GRM, which has already considered the allele frequency changes, hence
option `ong` is set to `true`.

See also [`randbrd`](@ref), [`aaocs`](@ref), [`iiocs`](@ref), [`iiocs`](@ref),
[`agocs`](@ref), [`igocs`](@ref).
"""
function ggocs(test, foo, bar, lmp, ngn, trait, fixed, plan, dF, F0; ε = 1e-6)
    @info "  - Directional selection GGOCS for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force = true)
    for ign = 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        G = grm(xy, lmp.chip, lmp.frq) + ε * I
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        g22 = G[ids, ids]
        ng = Select(ids, plan, ped, g22, trait, dF, ign; F0 = F0)
        reproduce!(ng, ped, xy, lmp, trait)
    end
    println()
    serialize("$test/$bar.ped", ped)
end

"""
    agocs(test, foo, bar, lmp, ngn, trait, fixed, noff, dF, F0; ε = 1e-6)
Optimal contribution selection with `G` relationship matrix for EBV, `A`
relationship matrix for constraint on `foo`.xy and `foo`.ped in directory `test`
for `ngn` generations. SNP linkage information are in DataFrame `lmp`. The
results are saved in `bar`.xy, `bar`.ped in directory `test`. The selection is
on a single trait `trait` with fixed effects `fixed`, which is a column name
vector in pedigree DataFrame. Parents are sampled according to `plan`. The
constraint ΔF is `dF`. `F0` is the inbreeding
coefficient of the `foo` population.

This function uses the TM1997 algorithm for OCS.

See also [`randbrd`](@ref), [`aaocs`](@ref), [`iiocs`](@ref), [`iiocs`](@ref),
[`ggocs`](@ref), [`igocs`](@ref).
"""
function agocs(test, foo, bar, lmp, ngn, trait, fixed, plan, dF, F0; ε = 1e-6)
    @info "  - Directional selection AGOCS for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force = true)
    for ign = 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        G = grm(xy, lmp.chip, lmp.frq) + ε * I
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        G = nrm(ped)
        g22 = G[ids, ids]
        ng = Select(ids, plan, ped, g22, trait, dF, ign; F0 = F0)
        reproduce!(ng, ped, xy, lmp, trait)
    end
    println()
    serialize("$test/$bar.ped", ped)
end

"""
    igocs(test, foo, bar, lmp, ngn, trait, fixed, plan, dF, F0; ε = 1e-6)
Optimal contribution selection with `G` relationship matrix for EBV, `IBD`
relationship matrix for constraint on `foo`.xy and `foo`.ped in directory `test`
for `ngn` generations. SNP linkage information are in DataFrame `lmp`. The
results are saved in `bar`.xy, `bar`.ped in directory `test`. The selection is
on a single trait `trait` with fixed effects `fixed`, which is a column name
vector in pedigree DataFrame. Parents are sampled according to `plan`. The
constraint ΔF is `dF`. `F0` is the inbreeding coefficient of the `foo`
population.

This function uses the TM1997 algorithm for OCS.

See also [`randbrd`](@ref), [`aaocs`](@ref), [`iiocs`](@ref), [`iiocs`](@ref),
[`ggocs`](@ref), [`agocs`](@ref).
"""
function igocs(test, foo, bar, lmp, ngn, trait, fixed, plan, dF, F0; ε = 1e-6)
    @info "  - Directional selection IGOCS for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force = true)
    for ign = 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        G = grm(xy, lmp.chip, lmp.frq) + ε * I
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        mid = size(ped, 1)
        g22 = irm(xy, lmp.chip, mid+1-length(ids):mid) + ε * I
        ng = Select(ids, plan, ped, g22, trait, dF, ign; F0 = F0)
        reproduce!(ng, ped, xy, lmp, trait)
    end
    println()
    serialize("$test/$bar.ped", ped)
end
