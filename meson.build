project(
    'gospl', 
    'c',
    version: '2024.09.01',
    meson_version: '>= 1.0',)

name = 'gospl'

add_languages('fortran')

platform = host_machine.system()
if platform == 'windows'
  add_project_link_arguments('-static', language: ['fortran', 'c'])
elif platform == 'darwin'
  add_project_link_arguments('-Wl,-rpath, "@loader_path"', language: ['fortran', 'c'])
else
  add_project_link_arguments('-Wl,-rpath,"$ORIGIN"', language: ['fortran', 'c'])
endif

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

install_subdir(name, install_dir: py.get_install_dir() / name, strip_directory: true)

fortran_source = custom_target('_fortranmodule.c',
  input : ['fortran/functions.pyf'],  
  output : ['_fortranmodule.c', '_fortran-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@']
)

py.extension_module('_fortran',
  ['fortran/functions.F90', fortran_source],
  incdir_f2py / 'fortranobject.c',
  subdir: 'gospl',
  include_directories: inc_np,
  dependencies: py_dep,
  install: true,
)