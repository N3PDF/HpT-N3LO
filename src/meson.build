sources = [
  'MellinFunc.cpp',
  'ThresExp.cpp',
  'ComplexDefs.cpp',
  'HarmonicSum.cpp',
  'SmallptExp.cpp'
]

higgslib = static_library(
  'higgspt',
  sources,
  include_directories : inc,
  dependencies : [
  	glib_dep,
	higgsfo_lib,
	cuba_lib,
	yaml_lib
  ],
  install: true
)
