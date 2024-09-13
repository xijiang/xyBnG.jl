"""
    extchr(base::AbstractString, dest::AbstractString, chr::Int...)
Extract chromosome `chr`s from directory `base` to `dest`.
"""
function extchr(base::AbstractString, dest::AbstractString, chr::Int...)
    base == dest && error("Can't write the results in the same directory of $base")
    isdir(dest) || mkpath(dest)
    name, nid = begin
        lines = readlines("$base/desc.txt")
        open("$dest/desc.txt", "w") do io
            println(io, lines[1])
            println(io, lines[2])
        end
        lines[1], parse(Int, lines[2])
    end
    lmp = deserialize("$base/$name.lmp")
    tmp = filter(row -> row.chr in chr, lmp)
    nlc = size(tmp, 1)
    hps = mapit("$base/$name.xy")
    open("$dest/$name.xy", "w") do io
        hdr = header("$base/$name.xy")
        write(io, Ref(hdr), [nlc, 2nid])
        write(io, hps[1:nlc, :])
    end
    serialize("$dest/$name.lmp", tmp)
end
