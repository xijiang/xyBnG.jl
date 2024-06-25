function tovcf(xy::AbstractMatrix, lmp::DataFrame, vcf::AbstractString)
    open(vcf, "w") do io
        println(io, "##fileformat=VCFv4.2")
        println(io, "##source=xy2vcf.jl")
        println(io, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
        println(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE")
        for (i, row) in enumerate(eachrow(xy))
            chrom = row[1]
            pos = row[2]
            ref = row[3]
            alt = row[4]
            println(io, "$chrom\t$pos\t.\t$ref\t$alt\t.\t.\t.\tGT\t0/1")
        end
    end
end
