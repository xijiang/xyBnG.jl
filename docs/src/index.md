```@meta
CurrentModule = xyBnG
```

# xyBnG

Documentation for [xyBnG](https://github.com/xijiang/xyBnG.jl).

```@index
```

```@autodocs
Modules = [xyBnG]
```

# File formats

## Linkage map

There are two kinds of maps. One is for the base population. This map has five
columns:

1. `chr::Int8`: Chromosome number.
2. `pos::Int32`: Mutation base pair position.
3. `ref::Char`: Ancestral base.
4. `alt::Char`: Mutation
5. `frq::Float64`: Allele frequency of 1 in the base population.

The other kind of map has 10+ columns, depends on number of traits:
1. `chr::Int8`: Chromosome number.
2. `pos::Int32`: Mutation base pair position.
3. `ref::Char`: Ancestral base.
4. `alt::Char`: Mutation
5. `frq::Float64`: Allele frequency of 1 in the founder population.
6. `chip::Bool`: Indicate if current SNP is on chips.
7. `dark::Bool`: Hidden and reference SNP.
8. `trait-1::Bool`: QTL for trait-1.
9. `trait-1_a::Float64`: Additive effects of trait-1 QTL
10. `trait-1_d::Float64`: Dominant effects of trait-1 QTL
11. Trait 2 info and on.

## `XY` files

An `XY` file stores a matrix and some extra information of the matrix. Each `XY`
file has 3 parts.
1. Header structure
    - `x::Int8` ≡ 'x'
    - `y::Int8` ≡ 'y'
    - `v::Int8` ≡ ' ', above are the magic chars for the `XY` file to be
      recognized.
    - `flus::Int8`, showing if the matrix is full, lower or upper triangle, or
      symmetric.
    - `major::Int8`, 0 for loci majored, which is default, and 1 for
      haplotype/ID majored.
    - `type::Int8`, element type of the matrix, defined in `_type`.
    - `r::Int8`, reserved
    - `u::Int8`,
        - 0 for SNP coding
        - 1 for IBD coding
        - 2 for genotype coding
        - 3 for `BitArray` coding
        - 4+ for else
2. Dimension, for the matrix stored.
3. Matrix part, written as binary.

# Package consistency check

There are a few functions to run to see if my program is working properly. 

## Package tests

Run `test` after the `xyBnG` environment is activated. It test

1. the availability of founder simulation programs, namely, `stdpopsim`,
   `tskit`, `slim`, `macs`
2. relationship matrix calculator, ARM, GRM, and IRM.
3. API for `xy` files.

## Functioning tests

1. `xyBnG.xps.a_tale_of_selection()`, tests
    - base simulation
    - founder sampling
    - random selection
    - directional selection
2. `xyBnG.xps.dosblup()`, compare all scenarios of selection schemes available
   with this package.