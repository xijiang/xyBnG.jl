"""
Provides functions to convert VCF files to and from `XY` files. All the VCF here
only deals with 2-allele SNP loci. For simulation purposes, there is no need to 
have complicated VCF files, at least at the moment.

2024-06-27
"""
module vcf

using DataFrames
using GZip
using Mmap
using Serialization
using Statistics
import xyBnG.Util: commas
using xyBnG.XY

function toxy(V::IO, X::AbstractString)
    @info "Counting valid loci and individuals in $V"
    nid = 0
    vld = Bool[]
    for line in eachline(V)
        if line[2] == '#'
            push!(vld, false)
            continue
        end
        if line[1] == '#' # split is expensive, so we do it only once
            nid = length(split(line)) - 9
            push!(vld, false)
            continue
        end
        m = 1
        for _ in 1:4
            m = findnext('\t', line, m) + 1
        end
        n = findnext('\t', line, m)
        push!(vld, n - m == 1)
    end
    
    nlc = sum(vld)
    hap = BitArray(undef, nlc, 2nid)
    lmp = DataFrame(
        chr = Int8[],
        pos = Int32[],
        ref = Char[],
        alt = Char[],
        frq = Float64[],
    )
    seekstart(V)
    oxy, omp = begin
        if length(X) > 3 && lowercase(X[end-2:end]) == ".xy"
            "$X", "$(X[1:end-3]).lmp"
        else
            "$X.xy", "$X.lmp"
        end
    end
    open(oxy, "w") do oo
        write(oo, Ref(XY.header(r = 1)), [nlc, 2nid])
        iln = ilc = 0
        for line in eachline(V)
            iln += 1
            vld[iln] || continue
            ilc += 1
            p = 1
            for _ in 1:5
                p = findnext('\t', line, p) + 1
            end
            chr, pos, _, ref, alt = split(line[1:p-1])
            push!(lmp, (parse(Int8, chr), parse(Int32, pos), ref[1], alt[1], 0.5))
            Alt = rand("01") # randomly flip the allele names
            #Alt = '1' # to test if read rightly.
            for j in 1:nid
                hap[ilc, 2j-1] = (line[p] == Alt)
                p += 2
                hap[ilc, 2j] = (line[p] == Alt)
                p += 2
            end
        end
        write(oo, hap)
        lmp.frq = vec(mean(hap, dims = 2))
        serialize(omp, lmp)
    end
end

function toxy(V::AbstractString, X::AbstractString)
    if length(V) > 4 && lowercase(V[end-3:end]) == ".vcf"
        open(V, "r") do ii
            toxy(ii, X)
        end
    elseif length(V) > 7 && lowercase(V[end-6:end]) == ".vcf.gz"
        GZip.open(V, "r") do ii
            toxy(ii, X)
        end
    else
        error("Unknown file type: $V")
    end
end

function toxy(V::AbstractString)
    p = findlast(".vcf", lowercase(V))
    if isnothing(p)
        error("Unknown file type: $V")
    end
    p = p[1] - 1
    toxy(V, V[1:p])
end

#=
"""
    dim(vcf::AbstractString)
Find the dimensions of a VCF file. Returns the number of loci and the number of
individuals.
"""
function dim(vcf::AbstractString)
    @info "Counting loci and individuals in $vcf"
    nlc, nid = 0, 0
    open(vcf, "r") do io
        for line in eachline(io)
            line[2] == '#' && continue
            if line[1] == '#'
                nid = length(split(line)) - 9
                continue
            end
            nlc += 1
            nlc % 100_000 == 0 && print("\r\t", commas("n_ID = $nid; n_Loci = $nlc"))
        end
    end
    println("\r\t", commas("n_ID = $nid; n_Loci = $nlc"))
    nlc, nid
end

"""
    firstline(vcf::AbstractString)
A temporary function to find the first line of a VCF file that is not a header.
"""
function firstline(vcf::AbstractString)
    open(vcf, "r") do io
        for line in eachline(io)
            line[1] ≠ '#' && return line
        end
    end
end

"""
    findnth(s::AbstractString, c::Union{AbstractString,AbstractChar}, n::Int)
Return the index of the first character of the nth occurrence of `c` in `s` and
the number of occurrences before the return. If the nth `c` is not found, return
`(nothing, i)`, where `i` is the number of accurances.
"""
function findnth(s::AbstractString, c::Union{AbstractString,AbstractChar}, n::Int)
    i, f = 0, 0
    while i < n
        t = findnext(c, s, f + 1)
        isnothing(t) && break
        f = t[1]
        i += 1
    end
    i < n && error("Only $i occurrences of $c found in $s")
    f, i
end

"""
    function line2v(line::AbstractString, av::Vector{Int8}; sep = '\t')
Convert a line of VCF file into a tuple of chromosome, position, ref and alt.
Also the frequency of allele 1. Write the alternative alleles into a vector of
Int8.
"""
function line2v(line::AbstractString, av::AbstractVector{Int8}; sep = '\t')
    f, _ = findnth(line, sep, 5)
    txt = split(line[1:f])
    chr = parse(Int8, txt[1])
    pos = parse(Int32, txt[2])
    ref = txt[4][1]
    alt = txt[5][1]
    f, _ = findnth(line, sep, 9)
    n, k = length(line), 0
    for i = f+1:4:n
        k += 1
        av[2k-1] = line[i]
        av[2k] = line[i+2]
    end
    av .-= 48  # Int('0') = 48

    return chr, pos, ref, alt, mean(av)
end

"""
    toxy(vcf::AbstractString, sxy::AbstractString; nln=10000)
Convert a VCF file into a `sxy.xy` and `sxy.lmp`. It deals with 10k loci
parallelly at a time. It also only deals with SNP VCF files. The SNPs should
have maximal 2 alleles only. Fields in the VCF file should be separated by tabs,
one and only one between two fields.
"""
function toxy(vcf::AbstractString, sxy::AbstractString; nln = 10000)
    nlc, nid = dim(vcf)
    hdr = XY.header(major = 1)
    lmp = DataFrame(
        chr = Int8[],
        pos = Int32[],
        ref = Char[],
        alt = Char[],
        frq = Float64[],
    )
    write(sxy * ".xy", Ref(hdr), [2nid, nlc])

    # Create blocks for parallel processing
    blks = collect(nln:nln:nlc)
    blks[end] < nlc && push!(blks, nlc)

    @info "Processing the genotypes:"
    open(sxy * ".xy", "w") do oo
        write(oo, Ref(hdr), [2nid, nlc])
        iloc, iblk, buf = 0, 1, IOBuffer()
        gt = zeros(Int8, 2nid, nln)
        for line in eachline(vcf)
            line[1] == '#' && continue
            ilc += 1
            println(buf, line)
            if iloc == blks[iblk]
                loci = String(take!(buf))
                iblk += 1
            end
        end     # skip header
    end
    return
    open(vcf, "r") do ii
        for line in eachline(ii)
            line[2] ≠ '#' && break
        end     # skip header
        ilc, jlc, ibk, buf = 0, 0, 1, String[]
        open(sxy * ".xy", "r+") do oo
            gt = Mmap.mmap(oo, Matrix{Int8}, (nlc, 2nid), 24)
            for line in eachline(ii)
                jlc += 1
                push!(buf, line)
                if jlc == blks[ibk]
                    Threads.@threads for i in eachindex(buf)
                        mmp[ilc+i, :] = vcfln2av(buf[i], view(gt, i + ilc, :))
                    end
                    print("\r\tProcessed loci: $(commas(jlc)) / $(commas(nlc))")
                    ibk += 1
                    ilc = jlc
                    empty!(buf)
                end
            end
        end
    end
    println('\n')
    serialize("$sxy.lmp", lmp)
end

function z2xy(zvcf::AbstractString, xy::AbstractString; nln = 10000)
    @info "Counting loci and individuals in $zvcf"
    nlc, nid = 0, 0
    open(`pigz -dc $zvcf`, "r+") do ii
        for line in eachline(ii)
            line[2] == '#' && continue
            if line[1] == '#'
                nid = length(split(line)) - 9
                continue
            end
            nlc += 1
            nlc % 100_000 == 0 &&
                print("\r\tn_ID = $(commas(nid)); n_Loci = $(commas(nlc))")
        end
        println("\r\tn_ID = $(commas(nid)); n_Loci = $(commas(nlc))\n")
    end

    open(`pigz -dc $zvcf`, "r+") do ii
        hdr = xyheader(nlc, 2nid)
        mmp = DataFrame(chr = zeros(Int8, nlc), pos = zeros(Int32, nlc), frq = zeros(nlc)) # map
        write(xy * "-hap.xy", Ref(hdr))

        bs = blksz(nlc, nln)
        blks = collect(bs:bs:nlc)
        blks[end] < nlc && push!(blks, nlc)

        @info "Processing the genotypes:"
        for line in eachline(ii)
            line[2] ≠ '#' && break
        end     # skip header
        ilc, jlc, ibk, buf = 0, 0, 1, String[]
        open(xy * "-hap.xy", "r+") do oo
            gt = Mmap.mmap(oo, Matrix{Int8}, (nlc, 2nid), 24)
            for line in eachline(ii)
                jlc += 1
                push!(buf, line)
                if jlc == blks[ibk]
                    Threads.@threads for i in eachindex(buf)
                        mmp[ilc+i, :] = vcfln2av(buf[i], view(gt, i + ilc, :))
                    end
                    print("\r\tProcessed loci: $(commas(jlc)) / $(commas(nlc))")
                    ibk += 1
                    ilc = jlc
                    empty!(buf)
                end
            end
        end
        println('\n')
        serialize("$xy-map.ser", mmp)
    end
end
=#

end # module VCF
