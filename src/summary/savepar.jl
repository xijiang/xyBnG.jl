"""
    savepar(dic, scenario)
Save the parameters of the current simulation into file `scenario`.
"""
function savepar(scenario, file)
    n = 0
    for k in keys(scenario)
        n = max(n, length(string(k)))
    end
    ios, iom = IOBuffer(), IOBuffer()

    for (k, v) in pairs(scenario)
        if typeof(v) ∈ (Cattle, Trait, Plan, Species)
            println(iom)
            println(iom, "[$k]")
            println(iom, v)
        elseif occursin("Tuple", string(typeof(v)))
            println(iom, "\n[Schemes]")
            for s in v
                println(iom, "         ", string(s))
            end
        else
            println(ios, lpad(k, n), ": ", v)
        end
    end
    open(file, "w") do io
        println(io, lpad("Started", n), ": ", time())
        # use Dates.unix2datetime(time()) to convert above time to DateTime
        print(io, String(take!(ios)))
        print(io, String(take!(iom)))
    end
    foo = joinpath(dirname(@__DIR__), "summary/summary.ipynb")
    bar = joinpath(dirname(file), "sums.ipynb")
    cp(foo, bar, force = true)
end
