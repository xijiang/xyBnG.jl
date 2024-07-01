"""
    phasedibd(env::AbstractString, vcf::AbstractString, lmp::AbstractString)

Calculate IBD relationships with genotypes in `vcf` and (plink formatted)
linkage map `lmp`, using Python package `phasedibd`. The trick here is that the
`phasedibd` package is old. It is using Cython language level 2 to compile. It
also needs older python environment. 

In my test, I used conda configured Python 3.7.16.

```bash
conda create -n py37 python=3.7
git clone https://github.com/23andMe/phasedibd
cd phasedibd
make
python setup.py install
```

You need to tell this Julia `phasedibd` function the environment name, in which
that Python `phasedibd` package was installed. For example, `py37`, as shown
above. Then the calculation is done by run ``Cmd`` `conda run -n py37 python
phasedibd.py`.`

"""
function phasedibd(env::AbstractString, vcf::AbstractString, lmp::AbstractString)
    @info "under construction"
end
