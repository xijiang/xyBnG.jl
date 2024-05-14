@testset "    Relationship matrices" begin
    rs = xyBnG.RS
    ped = DataFrame(id = 1:6, sire = [0, 0, 1, 1, 4, 5], dam=[0, 0, 2, 0, 3, 2])
    A = rs.nrm(ped)
    @test sum(A) == 18.625
    gt = [ 2  1  1  0  0  1  0  0  2  0  0  1  0  1
           0  0  1  0  1  1  0  1  0  0  1  0  0  0
           1  0  2  2  1  0  1  1  0  0  1  0  0  1
           1  0  1  1  2  1  1  0  0  1  0  0  1  1
           0  0  1  0  0  0  0  0  0  1  0  1  1  0
           0  2  0  1  0  2  2  1  1  2  1  1  2  2
           0  0  0  0  0  0  0  0  2  0  0  0  0  0
           2  2  2  2  2  2  2  2  2  2  2  2  2  1
           1  1  1  2  1  2  2  2  1  0  2  0  1  0
           2  0  2  1  2  1  0  0  2  0  1  0  0  0]
    G = rs.grm(gt)
    @test sum(G) < 1e-10
end
