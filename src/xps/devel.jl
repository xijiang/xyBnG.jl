function devel(; op = 1, nrpt = 1)
    if op == 1
        dosblup(nrpt = nrpt)
    elseif op == 2
        dir, fxy, fmp = "/mnt/a/store/xybng", "macs/macs.xy", "macs/macs.lmp"
        tdr = "cattle"
        milk = Trait("milk", 0.25, 10000)
        growth = Trait("growth", 0.5, 10000)
        sample_xy("$dir/$fxy", "$dir/$fmp", "$dir/$tdr", 200, 0.0, 50000, 10000, milk, growth)
    else
        @info "Good day!"
    end
end
