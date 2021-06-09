sources = [
  'MellinFunc.cpp',
  'ThresExp.cpp',
  'HighEnergyExp.cpp',
  'AnomalousDim.cpp',
  'ComplexDefs.cpp',
  'HarmonicSum.cpp',
  'SmallptExp.cpp',
  'CombinedRes.cpp'
]

higgslib = static_library(
  'higgspt',
  sources,
  include_directories : inc,
  dependencies : [
  	glib_dep,
	higgsfo_lib,
	cuba_lib,
	yaml_lib,
	cmpx_bessel
  ],
  install: true
)
