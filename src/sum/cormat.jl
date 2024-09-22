"""
    cormat(ma::AbstractMatrix, mb::AbstractMatrix)
Calculate the correlation between the off-diagonal elements of two matrices `ma`
and `mb`.
"""
function cormat(ma::AbstractMatrix, mb::AbstractMatrix)
    size(ma) == size(mb) || error("Matrix size not match")
    isequal(size(ma)...) || error("Not square matrices")
    va, vb = Float64[], Float64[]
    for i ∈ 1:size(ma, 1)-1
        append!(va, ma[i+1:end, i])
        append!(vb, mb[i+1:end, i])
    end
    cor(va, vb)
end

"""
    cormat(dir::AbstractString)
Calculate the correlation between the off-diagonal elements of the `A`, `G`, and
`IBD` matrices. G and IBD are calculated with the reference loci indicated in
the linkage map file.
"""
function cormat(dir::AbstractString)
    isfile("$dir/summary.ser") || error("No summary file found in $dir")
    ss = deserialize("$dir/summary.ser")
    cag, cai, cgi, npd = Float64[], Float64[], Float64[], ndigits(ss.repeat[end])
    for (irpt, cskm) ∈ eachrow(unique(select(ss, :repeat, :scheme)))
        @info "Processing repeat $irpt and scheme $cskm"
        tag = lpad(irpt, npd, '0')

        # with linkage map
        lmp = deserialize("$dir/$tag-founder.lmp")
        ref = lmp.dark
        frq = lmp.frq[ref]

        # with pedigree
        ped = deserialize("$dir/$tag-$cskm.ped")
        grt = repeat(ped.grt, inner=2)
        ugt = unique(grt)
        ngt = length(ugt)
        A = RS.nrm(ped)

        # with haplotypes
        hps = XY.mapit("$dir/$tag-$cskm.xy")
        for igt ∈ 1:ngt
            usnp = view(hps, ref, grt .== ugt[igt])
            D = RS.irm(usnp)
            gt = isodd.(usnp[:, 1:2:end]) + isodd.(usnp[:, 2:2:end])
            G = RS.grm(gt, p=frq)
            id = ped.grt .== ugt[igt]
            Aid = view(A, id, id)
            push!(cag, cormat(Aid, G))
            push!(cai, cormat(Aid, D))
            push!(cgi, cormat(G, D))
        end
    end
    ss.cag, ss.cai, ss.cgi = cag, cai, cgi
    serialize("$dir/summary.ser", ss)
end
