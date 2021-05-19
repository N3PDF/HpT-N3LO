#### Description

HpT-N3LO computes the approximation of the N3LO contribution to the Higgs
boson production from gluon fusion by combining the small-pt and threshold
resummations.


#### Dependency

`HpT-N3LO` relies on [HpT-MON](https://github.com/N3PDF/HpT-MON) for the computation
of the exact FO (up to NNLO) results. The instruction on how to install `HpT-MON` is
described in the `README` file.


#### Compile & run the code

Thanks to meson, compiling the code is straightforward:
```bash
meson setup builddir
cd builddir
meson compile
```

This will generate a binary called <kbd>higgs</kbd> in the `builddir` directory. To run
the code, use one of the run cards in the `runcards` folder as follows:
```bash
./higgs ../runcards/inputfile.yaml
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
