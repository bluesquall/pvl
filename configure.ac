# Copyright (C) 2011 M J Stanway <m.j.stanway@alum.mit.edu>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

AC_INIT([pvl],
        [0.0.1],
        [M J Stanway m.j.stanway@alum.mit.edu],
        [pvl])
AC_CONFIG_AUX_DIR([.build-aux])
AC_CONFIG_MACRO_DIRS([.build-aux/m4])
AM_INIT_AUTOMAKE([subdir-objects no-dist-gzip dist-xz -Wall -Werror foreign])

LF_CONFIGURE_CC
LF_CONFIGURE_CXX
AC_PROG_RANLIB
LF_HOST_TYPE
dnl LF_PROG_F77_PREFER_F2C_COMPATIBILITY
dnl LF_PROG_F77_PREFER_NATIVE_VERSION
dnl LF_PROG_F77
AC_PROG_FC
LF_SET_WARNINGS

AC_CONFIG_FILES([
   Makefile
   fortran/Makefile
])

AC_OUTPUT
