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

    nlc, nid, mid, rga, rgb = 1000, 10, 6, 1:6, 7:10
    loci, xy = ones(Bool, nlc), "mat.xy"
    M = rand(UInt8, nlc, 2nid)
    xyBnG.XY.mat(M, xy)
    hdr = xyBnG.XY.header(xy)
    hdr.u = 1
    xyBnG.XY.header!(xy, hdr)
    G = rs.irm(xy, loci, 1:nid)
    M = rs.irm(xy, loci, rga)
    M = rs.xirm(M, xy, loci, mid, nid)
    @test G == M
    rm(xy)
end
