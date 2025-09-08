project('pvl', 'fortran', meson_version: '>=0.49')

subdir('fortran')
executable('pvl',
  sources: pvl_sources,
  install: true)
install_data('fortran/example.inp',
  install_dir: get_option('datadir') / 'pvl')
