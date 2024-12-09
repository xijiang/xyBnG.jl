"""
    init!(xy::AbstractString, hdr::header, nrows::Int64, ncols::Int64)
Initialize an `xy` file with header `hdr` and dimensions `nrows` and `ncols`.
The matrix part is filled with zeros.
"""
function init!(xy::AbstractString, hdr::header, nrows::Int64, ncols::Int64)
    isfile(xy) && rm(xy, force = true)
    block = 24
    nbyte = sizeof(_type(hdr.type))
    if hdr.flus == Int8('F')
        block += nrows * ncols * nbyte
    else
        nrows == ncols || error("nrows must equal to ncols for full matrix")
        block += nrows * (nrows + 1) ÷ 2 * nbyte
    end

    if Sys.islinux() ||
       Sys.isbsd() ||
       Sys.isapple() ||
       Sys.isfreebsd() ||
       Sys.isnetbsd() ||
       Sys.isopenbsd() ||
       Sys.isunix()
        run(pipeline(`head -c $block /dev/zero`, xy))
    elseif Sys.iswindows()
        run(`fsutil file createnew $xy $block`)
    else
        error("Unsupported OS")
    end

    open(xy, "r+") do io
        write(io, Ref(hdr))
        write(io, [nrows, ncols])
    end
end
