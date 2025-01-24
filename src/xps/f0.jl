# Return the F0 with different relationship matrices. Before I used Falconer
# 1996, pp. 67. It is unfair for G, I and H matrices.

"""
    function meanoffd(a)
Possible the fastest way to calculate the mean off-diagonals of a square matrix.
"""
function meanoffd(a)
    m, n = size(a)
    m == n || error("Not square")
    (sum(a) - sum(diag(a))) / (n * (n - 1))
end

"""
    F0A(ped)
Empirical pedigree based inbreeding coefficient of the last generation of
pedigree `ped`.
"""
function F0A(ped)
    A = nrm(ped)
    #mean(diag(A)[ped.grt .== ped.grt[end]]) - 1
    lg = ped.grt .== ped.grt[end]
    T = view(A, lg, lg)
    meanoffd(T)
end

"""
    F0H2(xy, lmp, ped)
Empirical inbreeding coefficient of the last generation of pedigree `ped` based
on the heterozygosity differences between the first and last generation of chip
loci.
"""
function F0H2(xy, lmp, ped)
    hap = Int8.(isodd.(XY.mapit(xy)))
    fg = repeat(ped.grt .== ped.grt[begin], inner = 2)
    f0 = vec(mean(hap[lmp.chip, fg], dims = 2))
    h0 = Sum.v_snphet(f0)
    vld = h0 .> 0.
    lg = repeat(ped.grt .== ped.grt[end], inner = 2)
    ft = vec(mean(hap[lmp.chip, lg], dims = 2))
    ht = Sum.v_snphet(ft)
    mean((h0[vld] .- ht[vld]) ./ h0[vld])
end

"""
    F0H(xy, lmp, ped)
Empirical inbreeding coefficient of the last generation of pedigree `ped` based
on the genomic relationship matrix with q0 .= 0.5 and with chip loci.
"""
function F0H(xy, lmp, ped)
    hap = Int8.(isodd.(XY.mapit(xy)))
    id = ped.id[ped.grt .== ped.grt[end]]
    gt = hap[lmp.chip, 2id .- 1] + hap[lmp.chip, 2id]
    H = grm(gt, p = ones(size(gt, 1)) * 0.5)
    meanoffd(H)
end

"""
    F0G(xy, lmp, ped)
Empirical inbreeding coefficient of the last generation of pedigree `ped` based
on the genomic relationship matrix with chip loci.
"""
function F0G(xy, lmp, ped)
    hap = Int8.(isodd.(XY.mapit(xy)))
    gt = hap[lmp.chip, 1:2:end] + hap[lmp.chip, 2:2:end]
    frq = mean(gt[:, ped.grt .== ped.grt[begin]], dims = 2) / 2
    G = grm(gt; p = vec(frq))
    meanoffd(G)
end

"""
    F0I(xy, lmp, ped)
Empirical inbreeding coefficient of the last generation of pedigree `ped` based
on the identity by descent matrix with chip loci.
"""
function F0I(xy, lmp, ped)
    lg = repeat(ped.grt .== ped.grt[end], inner = 2)
    hap = XY.mapit(xy)
    lhp = view(hap, lmp.chip, lg) #hap[lmp.chip, lg]
    T = irm(lhp)
    meanoffd(T)
end
