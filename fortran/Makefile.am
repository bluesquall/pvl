# Copyright (C) 1988 Eleftherios Gkioulekas <lf@amath.washington.edu>
#  
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 
SUFFIXES = .f90
.f.o:
	$(F77) $(FFLAGS) -c $<

bin_PROGRAMS = pvl
pvl_SOURCES = splutil.f90 Wrench.f90 Simeqn.f90 Pvlmod.f90 Forces.f90 Pvl.f90

pkgdata_DATA = example.inp
