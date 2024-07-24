using xyBnG
using Documenter

DocMeta.setdocmeta!(xyBnG, :DocTestSetup, :(using xyBnG); recursive = true)

makedocs(;
    modules = [xyBnG],
    authors = "Xijiang Yu <xijiang@users.noreply.github.com> and contributors",
    sitename = "xyBnG.jl",
    format = Documenter.HTML(;
        canonical = "https://xijiang.github.io/xyBnG.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/xijiang/xyBnG.jl", devbranch = "main")
