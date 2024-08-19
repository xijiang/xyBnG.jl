function devel()
    base = "rst/tskit"
    test = "rst/devel"
    isdir(test) || mkpath(test)
    sample_xy(
        "$base/BosTau.xy",
        "$base/BosTau.lmp",
        test,
        200,
        0.0,
        3000,
        1000,
        Trait("milk", 0.25, 1000),
    )
end
