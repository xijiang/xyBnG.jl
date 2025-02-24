"""
    resum(dir::AbstractString)
This is an amendment to the `xysum` function in directory `dir`, where the data
are kept. This is needed when the formulae of some indices are changed.
"""
function resum(dir::AbstractString, trait::Trait)
    isfile("$dir/summary.ser") || error("No summary file found in $dir")
    ss = deserialize("$dir/summary.ser")
    nrpt = ss.repeat[end]
    mv("$dir/summary.ser", "$dir/summary.bak", force = true)
    for (irpt, cskm) ∈ eachrow(unique(select(ss, :repeat, :scheme)))
        @info "Processing repeat $irpt and scheme $cskm"
        tag = lpad(irpt, ndigits(nrpt), '0')
        ped = deserialize("$dir/$tag-$cskm.ped")
        lmp = deserialize("$dir/$tag-founder.lmp")
        xy = "$dir/$tag-$cskm.xy"
        st = xysum(ped, xy, lmp, trait)
        savesum("$dir/summary.ser", st)
    end
end
