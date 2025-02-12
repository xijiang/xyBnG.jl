module xy
using DataFrames
using Serialization
using xyBnG.XY

function toplink() end

function togatc() end

"""
    tovcf(fxy::AbstractString, fmp::AbstractString, vcf::AbstractString)
Convert an `XY` file to a VCF file.
This function will only output the chip SNP loci.
"""
function tovcf(fxy::AbstractString, fmp::AbstractString, vcf::AbstractString)
    lmp = deserialize(fmp)
    nlc, nhp = XY.dim(fxy)
    nid = nhp ÷ 2
    hdr = XY.header(fxy)
    gt = hdr.u == 1 ? isodd.(XY.mapit(fxy)) : XY.mapit(fxy)
    open(vcf, "w") do io
        println(io, "##fileformat=VCFv4.2")
        println(io, "##source=xyBnG.Conn.xy.toxy")
        println(io, "##FILTER=<ID=PASS,Description=\"All filters passed\">")
        println(io, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
        print(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for i = 1:nid
            print(io, "\txy-" * lpad(string(i), ndigits(nid), '0'))
        end
        println(io)
        for i = 1:nlc
            lmp.chip[i] || continue
            buf = IOBuffer()
            print(
                buf,
                lmp.chr[i],
                '\t',
                lmp.pos[i],
                '\t',
                i,
                '\t',
                lmp.ref[i],
                '\t',
                lmp.alt[i],
                "\t.\tPASS\t.\tGT",
            )
            for j = 1:nid
                print(buf, '\t', join(Int8.(gt[i, 2j-1:2j]), '|'))
            end
            println(io, String(take!(buf)))
        end
    end
end

"""
    togt(
        fxy::AbstractString,
        fmp::AbstractString,
        fgt::AbstractString,
        id; # id::Vector{Int}
        snp = :chip,
    )
Convert the `id` and SNP class `snp` in an `XY` file to a Zooroh file in `GT`
format, which has 5+ columns:

1. chromosome name, named as `chr1`, `chr2`, etc.
2. marker position, in bp
3. reference allele
4. alternative allele
5. genotype of the first individual, and so on.

Argument `snp` can be one of the following:

- `:chip`: only output the chip SNPs
- `:dark`: only output the reference SNPs
- `:trait`: e.g., `:growth`, output the QTL for trait `growth`
"""
function togt(
    fxy::AbstractString,
    fmp::AbstractString,
    fgt::AbstractString,
    id; # id::Vector{Int}
    snp = :chip,
)
    nlc, nhp = XY.dim(fxy)
    lmp = deserialize(fmp)
    nlc == nrow(lmp) || throw(ArgumentError("The number of loci in the XY file is not equal to the number of loci in the map file."))
    cls = names(lmp)
    hasproperty(lmp, snp) || throw(ArgumentError("The SNP class `$snp` is not found in the map file."))
    olc = lmp[!, snp]

    ph = 2id .- 1  # paternal haplotypes
    mh = 2id       # maternal haplotypes

    hps = XY.mapit(fxy)
    hdr = XY.header(fxy)
    gt = if hdr.u == 1
        isodd.(view(hps, olc, ph)) + isodd.(view(hps, olc, mh))
    else
        view(hps, olc, ph) + view(hps, olc, mh)
    end
    open(fgt, "w") do io
        df = view(lmp, olc, [:chr, :pos, :ref, :alt])
        for i in 1:nrow(df)
            print(io, "chr", join(df[i, :], ' '), ' ')
            println(io, join(gt[i, :], ' '))
        end
    end
end

end # module xy
