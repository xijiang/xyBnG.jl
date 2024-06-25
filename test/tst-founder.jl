@testset "   Requirements 4 fdr sim" begin
    # These python package and software are required for the founder simulation
    @test !isnothing(Sys.which("stdpopsim"))
    @test !isnothing(Sys.which("tskit"))
    @test !isnothing(Sys.which("slim"))
    @test !isnothing(Sys.which("macs"))
end
