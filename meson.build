
project('coloquinte', 'cpp',
  version : '0.3.1',
  license : 'MIT',
  default_options: [
     'cpp_std=c++17',
  ],
)

coloquinte_includes = include_directories('src')
boost_dep = dependency('boost', required: true)
lemon_dep = dependency('lemon', required: true)
eigen_dep = dependency('eigen3', required: true)
thread_dep = dependency('threads', required: true)

sources = [
  'src/coloquinte.cpp',
  'src/parameters.cpp',
  'src/place_global/net_model.cpp',
  'src/place_global/density_legalizer.cpp',
  'src/place_global/density_grid.cpp',
  'src/place_global/place_global.cpp',
  'src/place_detailed/legalizer.cpp',
  'src/place_detailed/abacus_legalizer.cpp',
  'src/place_detailed/tetris_legalizer.cpp',
  'src/place_detailed/row_legalizer.cpp',
  'src/place_detailed/place_detailed.cpp',
  'src/place_detailed/detailed_placement.cpp',
  'src/place_detailed/incr_net_model.cpp',
  'src/place_detailed/row_neighbourhood.cpp',
  'src/place_global/transportation.cpp',
  'src/place_global/transportation_1d.cpp',
]

libcoloquinte = shared_library('coloquinte',
                               sources,
                               include_directories : coloquinte_includes,
                               dependencies : [boost_dep, lemon_dep, eigen_dep, thread_dep],
                               install : true)

install_headers('src/coloquinte.hpp')

coloquinte_dep = declare_dependency(link_with : libcoloquinte, include_directories: coloquinte_includes)

