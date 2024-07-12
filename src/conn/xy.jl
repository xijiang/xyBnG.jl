module xy
using DataFrames
using Serialization
using xyBnG.XY

function toplink()
end

function togatc()
end

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
    gt = hdr.u == 1 ?
        isodd.(XY.mapit(fxy)) :
        XY.mapit(fxy)
    open(vcf, "w") do io
        println(io, "##fileformat=VCFv4.2")
        println(io, "##source=xyBnG.Conn.xy.toxy")
        println(io, "##FILTER=<ID=PASS,Description=\"All filters passed\">")
        println(io, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
        print(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for i in 1:nid
            print(io, "\txy-" * lpad(string(i), ndigits(nid), '0'))
        end
        println(io)
        for i in 1:nlc
            lmp.chip[i] || continue
            buf = IOBuffer()
            print(buf, lmp.chr[i], '\t', lmp.pos[i], '\t', i, '\t', lmp.ref[i], '\t', lmp.alt[i], "\t.\tPASS\t.\tGT")
            for j in 1:nid
                print(buf, '\t', join(Int8.(gt[i, 2j-1:2j]), '|'))
            end
            println(io, String(take!(buf)))
        end
    end
end

end # module xy