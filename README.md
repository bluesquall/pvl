# PVL - Propeller Vortex Lattice

This repository provides an example implementation of the Vortex Lattice
Method (VLM) for marine propeller design. The algorithms are public domain.

## Install
This project uses [`meson`][meson] to build. For example:
```
meson setup /tmp/build/pvl
meson compile -C /tmp/build/pvl
sudo meson install -C /tmp/build/pvl
```

## Quick-start
The installation instructions above will place an example PVL input file
`example.inp` in `/usr/local/share/pvl`. To use that input file, start `pvl`
and provide the path to that input file when prompted. Alternatively:
```
echo "/usr/local/share/pvl/example.inp" | pvl
```

## Source
The Fortran source is based on J. E. Kerwin's PVL implementation in the
lecture notes for the MIT course 13.04/2.23 Hydrofoils and Propellers.
[MIT OpenCourseWare (OCW)][ocw] hosts subject material for the course as
taught [Fall 2003][dspace-13.04-2003] and [Spring 2007][ocw-2.23-2007].

Example implementations in more languages are in the works.

This repository is intended as an educational resource.

## License
Source is licensed under the [GPLv2][GPLv2] -- see `LICENSE`.

## see also
For a more advanced propeller design tool (based on the same basic
principles) refer to [the OpenProp project][openprop].

_____________
_____________
[ocw]: https://ocw.mit.edu 
[dspace-13.04-2003]: https://dspace.mit.edu/bitstream/handle/1721.1/36898/13-04Fall2003/OcwWeb/Ocean-Engineering/13-04Fall2003/CourseHome/index.htm?sequence=1
[ocw-2.23-2007]: https://ocw.mit.edu/courses/mechanical-engineering/2-23-hydrofoils-and-propellers-spring-2007/
[openprop]:  https://epps.com/openprop
[GPLv2]: https://spdx.org/licenses/GPL-2.0.html
[meson]: https://mesonbuild.com
