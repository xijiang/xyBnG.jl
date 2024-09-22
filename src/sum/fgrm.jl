"""
    fgrm(dir::AbstractString)
Update `summary.ser` with `F_grm`, which are calculated from the `xy` files in
directory `dir`. GRM is calculated with the chip loci.
"""
function fgrm(dir::AbstractString)
    isfile("$dir/summary.ser") || error("No summary file found in $dir")
    ss = deserialize("$dir/summary.ser")
    F, npd = Float64[], ndigits(ss.repeat[end])
    for (irpt, cskm) ∈ eachrow(unique(select(ss, :repeat, :scheme)))
        @info "Processing repeat $irpt and scheme $cskm"
        tag = lpad(irpt, npd, '0')

        # with linkage map
        lmp = deserialize("$dir/$tag-founder.lmp")
        chp = lmp.chip
        frq = lmp.frq[chp]

        # with pedigree
        ped = deserialize("$dir/$tag-$cskm.ped")
        grt = repeat(ped.grt, inner=2)
        ugt = unique(grt)
        ngt = length(ugt)
        tf = zeros(ngt)

        # with haplotypes
        hps = isodd.(XY.mapit("$dir/$tag-$cskm.xy"))
        for igt ∈ 1:ngt
            snps = view(hps, chp, grt .== ugt[igt])
            gt = snps[:, 1:2:end] .+ snps[:, 2:2:end]
            G = RS.grm(gt, p=frq)
            tf[igt] = mean(diag(G) .- 1)
        end
        append!(F, tf)
    end
    ss.fgrm = F
    serialize("$dir/summary.ser", ss)
end
