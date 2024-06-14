@testset "          XY module tests" begin
    t0, t1 = Int8(0), Int8(1)
    xy = xyBnG.XY
    # Test _type function
    @test xy._type(Int8) == 1
    @test xy._type(t1) == Int8
    
    # Test header struct
    hdr = xy.header()
    hdr.x, hdr.y, hdr.v, hdr.flus == Int8('x'), Int8('y'), Int8(' '), Int8('F')
    hdr.major, hdr.type, hdr.r, hdr.u == t0, t1, t0, t0
    
    hdr2 = xy.header(flus = 'F', major = 0, type = 1, u = 0)
    @test xy.isEqual(hdr, hdr2)

    # Test init! function
    xy_file = "test.xy"
    hdr = xy.header()
    nrows, ncols = 10, 5
    xy.init!(xy_file, hdr, nrows, ncols)
    @test isfile(xy_file)
    @test xy.isEqual(xy.header(xy_file), hdr)
    @test xy.dim(xy_file) == (nrows, ncols)
    rm(xy_file)

    # Test sub function
    ixy_file, oxy_file = "input.xy", "output.xy"
    xy.init!(ixy_file, hdr, nrows, ncols)
    rows, cols = [1, 3, 5], [2, 4]
    xy.sub(ixy_file, rows, cols, oxy_file)
    @test isfile(oxy_file)
    @test xy.dim(oxy_file) == (length(rows), length(cols))
    rm.([ixy_file, oxy_file])
    
    # Test mat function
    m = [1 2 3; 4 5 6; 7 8 9]
    mat_file = "matrix.xy"
    xy.mat(m, mat_file)
    @test isfile(mat_file)
    @test xy.dim(mat_file) == size(m)
    @test xy.mat(mat_file) == m
    rm(mat_file)
    
    # Test append! function
    xy.init!(ixy_file, hdr, nrows, ncols)
    xy.init!(oxy_file, hdr, nrows, ncols)
    xy.append!(ixy_file, oxy_file)
    @test xy.dim(ixy_file) == (nrows, 2 * ncols)
    rm.([ixy_file, oxy_file])
    
    # Test merge function
    jxy_file = "input2.xy"
    xy.init!(ixy_file, hdr, nrows, ncols)
    xy.init!(jxy_file, hdr, nrows, ncols)
    xy.merge(ixy_file, jxy_file, oxy_file)
    @test isfile(oxy_file)
    @test xy.dim(oxy_file) == (nrows, 2ncols)
    xy.merge(ixy_file, jxy_file, oxy_file, horizontal = false)
    @test xy.dim(oxy_file) == (2nrows, ncols)
    rm.([ixy_file, jxy_file, oxy_file])

    # Test code function
    nrows, ncols = 10, 6
    m1 = rand(Int8.([0, 1]), nrows, ncols)
    write(ixy_file, Ref(hdr), [nrows, ncols], m1)
    xy.code(ixy_file, oxy_file)
    m2 = xy.mat(oxy_file)
    coded = true
    for i in eachrow(m)
        length(unique(i)) == ncols && coded
    end
    @test coded
    @test m1 == isodd.(m2)
    rm.([ixy_file, oxy_file])
    
    # Test transpose! function
    write(ixy_file, Ref(hdr), [nrows, ncols], m1)
    xy.transpose!(ixy_file, oxy_file)
    m2 = xy.mat(oxy_file)
    @test m1 == m2'
    @test xy.dim(oxy_file) == (ncols, nrows)
    rm.([ixy_file, oxy_file])
    
    # Test snp2gt function
    write(ixy_file, Ref(hdr), [nrows, ncols], m1)
    xy.snp2gt(ixy_file, oxy_file)
    m2 = xy.mat(oxy_file)
    @test m1[:, 1:2:end] + m1[:, 2:2:end] == m2
    rm.([ixy_file, oxy_file])
    
    @info "Under construction: transpose! for half-matrix, and bed function"
end
