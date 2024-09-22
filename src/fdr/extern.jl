"""
    ts_base(pop::Cattle, dir::AbstractString)

Simulate a cattle population of name `pop.name`, and `pop.nid` ID in `dir`.

Note:
- This is a simulation with the coancestor/backward simulator msprime.
- Needs to have `tskit`, `msprime`, `scipy` and `stdpopsim` installed.
"""
function ts_base(pop::Cattle, dir::AbstractString)
    '~' ∈ dir && error("Character '~' is forbidden in dir")
    isdir(dir) || mkpath(dir)
    @info "Simulating a cattle population of name $(pop.name), and $(pop.nid) ID in $dir
      Note: This is a simulation with the coancestor/backward simulator msprime."
    Threads.@threads for chr = 1:29
        print(" $chr")
        cmd = pipeline(
            `stdpopsim BosTau -c $chr -o $dir/$chr.ts
 -d HolsteinFriesian_1M13 Holstein_Friesian:$(pop.nid)`,
            stderr=devnull,
        )
        run(cmd)
    end
    open("$dir/desc.txt", "w") do io
        println(io, pop.name)
        println(io, pop.nid)
    end
end

"""
    ts_base(pop::Species, dir::AbstractString)

Default function to catch unsupported species.
"""
function ts_base(pop::Species, dir::AbstractString)
    @info "This species is not supported yet." # or as below
    # https://popsim-consortium.github.io/stdpopsim-docs/stable/api.html#sec-api-generic-models
end

"""
    macs_base(pop::Cattle, dir::AbstractString)
Simulate a populatoin with `MaCS` into `dir`. The parameters are adapted from
``https://academic.oup.com/g3journal/article/2/4/425/6026056?login=true#supplementary-data``.
This simulation used the `Ne = 100` one.

Note, the command `macs` needs to be in a searchable path with no space
character in it. Files like `chr.1`, `log.1` are generated.
"""
function macs_base(pop::Cattle, dir::AbstractString)
    nid = pop.nid
    isdir(dir) || mkpath(dir)
    Ne = 100
    μ = 2.5e-8 * (4Ne)
    # [Chromosome length in bp](https://www.ncbi.nlm.nih.gov/projects/r_gencoll/ftp_service/nph-gc-ftp-service.cgi/?HistoryId=MCID_642d5d40ceff2e2c64293c60&QueryKey=1&ReleaseType=RefSeq&FileType=GENOME_FASTA&Flat=true)
    chr = [
        158534110,
        136231102,
        121005158,
        120000601,
        120089316,
        117806340,
        110682743,
        113319770,
        105454467,
        103308737,
        106982474,
        87216183,
        83472345,
        82403003,
        85007780,
        81013979,
        73167244,
        65820629,
        63449741,
        71974595,
        69862954,
        60773035,
        52498615,
        62317253,
        42350435,
        51992305,
        45612108,
        45940150,
        51098607,
    ]
    r = 1e-8 * 4Ne
    # Scaled time
    tm =
        [
            10,
            25,
            50,
            100,
            200,
            300,
            400,
            500,
            600,
            700,
            800,
            900,
            1000,
            2000,
            3000,
            4000,
            5000,
            6000,
            7000,
            8000,
            9000,
            10000,
            20000,
            40000,
            60000,
            80000,
            100000,
            200000,
            400000,
            600000,
            800000,
        ] / (4Ne)
    # Scaled population size
    ps =
        [
            175,
            200,
            350,
            500,
            700,
            820,
            850,
            900,
            1000,
            1100,
            1275,
            1300,
            1200,
            2000,
            2500,
            3000,
            3200,
            3500,
            3800,
            4000,
            4200,
            4500,
            5456,
            7367,
            9278,
            11190,
            13101,
            22658,
            41772,
            60886,
            80000,
        ] / Ne
    eN = ""
    for i ∈ eachindex(tm) # = 1:length(tm)
        eN *= " -eN $(tm[i]) $(ps[i])"
    end
    macs = Sys.which("macs")
    (isnothing(macs) || any(isspace.(collect(macs)))) && error("Command `macs` error")
    @info "  - Simulating a cattle population with MaCS into $dir"
    Threads.@threads for i ∈ eachindex(chr) # = 1:length(chr)
        print(" $i")
        cmd = "$macs $(2nid) $(chr[i]) -t $μ -r $r" * eN
        cmd = Cmd(convert(Vector{String}, split(cmd)))
        run(pipeline(cmd, stdout="$dir/chr.$i", stderr="$dir/log.$i"))
    end
    open("$dir/desc.txt", "w") do io
        println(io, pop.name)
        println(io, pop.nid)
    end
    println()
end
