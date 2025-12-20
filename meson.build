
project('coloquinte', 'cpp',
  version : '0.3.1',
  license : 'MIT',
  default_options: [
     'cpp_std=c++17',
  ],
)

if build_machine.system() == 'darwin'
  add_project_arguments('-mmacosx-version-min=13.0', language: ['c','cpp'])
endif

coloquinte_includes = include_directories('src')
boost_base_dep = dependency('boost', required: true)
if boost_base_dep.version().split('.')[1].to_int() < 89
  boost_modules = ['system', 'filesystem', 'iostreams', 'program_options', 'unit_test_framework']
else
  boost_modules = ['filesystem', 'iostreams', 'program_options', 'unit_test_framework']
endif
boost_dep = dependency('boost', required: true, modules: boost_modules)

eigen_dep = dependency('eigen3', required: true)
thread_dep = dependency('threads', required: true)

cmake =import('cmake')
cmake_opts = cmake.subproject_options()
cmake_opts.add_cmake_defines({'CMAKE_POSITION_INDEPENDENT_CODE': true})
cmake_opts.set_install(false)
lemon = cmake.subproject('lemon', options: cmake_opts)
lemon_dep = lemon.dependency('lemon')

sources = [
  'src/coloquinte.cpp',
  'src/parameters.cpp',
  'src/export.cpp',
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

tests = [
    'test_net_model'
,   'test_density_legalizer'
,   'test_density_grid'
,   'test_row_legalizer'
,   'test_legalizer'
,   'test_incr_net_model'
,   'test_detailed_placement'
,   'test_row_neighbourhood'
,   'test_transportation'
,   'test_expansion'
]

foreach testcase: tests
  test_exec = executable(testcase, 'test/' + testcase + '.cpp', dependencies: [boost_dep, coloquinte_dep])
  test(testcase, test_exec)
endforeach
