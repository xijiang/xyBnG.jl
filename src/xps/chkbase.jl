"""
    chkbase(dir::AbstractString, pop::Species; base = ts_base)
Check if the base population `pop` is already created. If not, create it.
"""
function chkbase(dir::AbstractString, pop::Species; base = ts_base)
    if isfile("$dir/desc.txt")
        desc = readlines("$dir/desc.txt")
        desc[1] == pop.name && parse(Int, desc[2]) ≥ pop.nid && return
    end
    if base == ts_base
        ts_base(pop, dir)
        Conn.TS.toxy(dir)
    elseif base == macs_base
        macs_base(pop, dir)
        Conn.MaCS.toxy(dir)
    else
        error("Unknown base population creator: $base")
    end
end
