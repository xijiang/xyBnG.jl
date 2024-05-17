function ocs_sel_milk(; bar = "ocs")
    dir, fixed = "rst/tale", ["grt"]

    for ext in ["ped", "xy"]
        cp("$dir/rand.$ext", "$dir/$bar.$ext", force=true)
    end

    @info "Directional selection test: index selection"
    milk   = Trait("milk",   0.3, 10000; sex = 0, vd = 0.1)
    growth = Trait("growth", 0.5, 10000; sex = 2, vd = 0.2)
    lmp = deserialize("$dir/founder.lmp")
    ped = deserialize("$dir/$bar.ped")
    xy  = "$dir/$bar.xy"
    plan = Plan(25, 50, 200, mate = :random) # default :hierarchical

    for igrt in 1:10 # AABLUP
        println()
        @info lpad("<--- Generation $igrt/10 --->", 40)
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, milk, growth)
        giv, R = begin
            G = nrm(ped)
            giv = inv(G)
            R = G[ids, ids]
            giv, R
        end
        predict!(ids, ped, fixed, giv, milk)
        tgt = select(ped, "ebv_milk" => "idx", "sex")[ids, :]
        K = 0.09
        c = TM1997(tgt, R, K)
        #= with Theo's or
        dat = Float64.([tgt.idx  tgt.sex .== 1 tgt.sex .== 0])
        cb = fungencont(dat, R, K)
        return [c cb]  # passed
        =#
        ng = Select(ids, plan, ped, c)
        reproduce!(ng, ped, xy, lmp, milk, growth)
    end
    serialize("$dir/$bar.ped", ped)
end

function sumtest()
    n, s, x = 200, 0, 0
    for i in 1:1000
        t = rand(n)
        t /= sum(t)
        y = sum(Int.(round.(t * n)))
        y > n && (x += 1)
        y < n && (s += 1)
    end
    println(s, " ", x)
end

## ToDo
# 1. unify Select
# 2. finish this function
# 3. Dos ocs.
# 4. More simulation on dos ocs
# 5. Oda's scheme.