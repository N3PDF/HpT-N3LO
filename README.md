[![DOI](https://zenodo.org/badge/368804644.svg)](https://zenodo.org/badge/latestdoi/368804644)

#### Descriptioh

HpT-N3LO computes the approximation of the N3LO contribution to the Higgs
boson production from gluon fusion by combining the small-pt and threshold
resummations.


#### Dependency

`HpT-N3LO` relies on the following libraries:

* [HpT-MON](https://github.com/N3PDF/HpT-MON): for the computation of the exact 
  FO (up to NNLO) results. The instruction on how to install `HpT-MON` is described 
  in the `README` file.

* [Complex Bessel](https://blog.joey-dumont.ca/complex_bessel/): for the computation
  of various Bessel functions with complex arguments. This is not currently being
  used yet but might be useful in the future (plus stripping it away from the rest
  of the code also demands some time).


#### Compile & run the code

Thanks to meson, compiling the code is straightforward:
```bash
meson setup builddir
cd builddir
meson compile
```

This will generate two executables called `higgs-pt` and `higgs-n` in the `builddir` 
directory. To run the code, use one of the run cards in the `runcards` folder as follows:
```bash
./higgs-n ../runcards/Mellin-HpT-as-N.yaml    (for results as a function of N)
./higgs-pt ../runcards/Mellin-HpT-as-pt.yaml  (for results as a function of pt)
```

Every time changes are made, the code can be re-compiled by running `meson compile`
inside the `builddir` directory.

Finally, in case one wants to install the header files and library system-wide, this
can be done by running the following:
```bash
meson install
```
This, by default, will install the header files in `/{prefix}/include/higgs-fo` and
add `higgsfo.pc` to the PKG path.
