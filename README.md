OpenMalaria
============

This git-repository contains source-code for OpenMalaria, a simulator program
for studying malaria epidemiology and the impacts of interventions against
malaria.

------

For further documentation, take a look at our
[wiki](https://github.com/SwissTPH/openmalaria/wiki).

Also find our code on zenodo, with DOI-references:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10534022.svg)](https://doi.org/10.5281/zenodo.10534022)



Schema documentation for the XML input files can be found
[on this website](https://swisstph.github.io/openmalaria/)
or in the `schema` folder (both past releases and development version).

You can download the latest build here:
[releases](https://github.com/SwissTPH/openmalaria/releases)

Version
======

The schema version is specified in the following places (all need to be updated
when releasing a new version):

*   DocumentLoader.h
*   schema/scenario.xsd, demography.xsd, etc. (all XSD files without a version number)
*   schema/CMakeLists.txt (namespace-map)
*   copy build/schema/scenario_current.xsd to schema/scenario_XX.xsd
*   test/*.xml — update http://openmalaria.org/schema/scenario_XX and (optionally) schemaVersion="XX"
*   version.txt — needed for build service

In theory the "schema namespace version" doesn't need to match the "OpenMalaria"
version and we could update the latter without requiring changes to XML files,
however currently we keep both synchronised (in some ways this is simpler).


Installation instructions:
==================
See [INSTALLATION](https://github.com/SwissTPH/openmalaria/wiki/UserGuide).

Build instructions:
===================

```
mkdir build && cd build
ccmake ..
Press 'c', look over options, press 'c' again and 'g'
make -j4
ctest -j4
```

For testing and development, ideally use debug builds (which enable some
asserts to do with simulation time usage).

Code subdirectories:
=============
|- dir    -|- description -|
|----------|:-------------------------------------------------------------------------:|
| contrib | Third-party libraries, distributed under the same repo for convenience.   |
| model   | Source code for the malaria model.                                        |
| test    | High-level testing: test scenarios with expected outputs. Also run-time files: densities.csv, scenario_?.xsd, Nv0scenario*.txt. |
| unittest| Low-level testing: unittests for the model using cxxunit. |
| util    | Extra scripts associated with OpenMalaria. |
| schema  | scenario schema files (see schema/policy.txt for details) |
| schema/scenario.xsd | The latest (partial) schema file. |
| schema/entomology.xsd, schema/demography.xsd, etc | components of the latest schema, included from scenario.xsd. |
| schema/scenario_*.xsd | Copies of released schema versions, with all components inlined in the same file. |

This git repository is currently maintained by members of the [Global Disease Modelling Team](https://www.thekids.org.au/our-research/early-environment/infection-and-vaccines/global-disease-modelling/) of __The Kids Research Institute Australia__,  the [Disease Modelling Group](https://www.swisstph.ch/en/about/eph/disease-modelling/) of the __Swiss Tropical and Public Health institute__ along with the [SciCORE Center for Scientific Computing, University of Basel](https://scicore.unibas.ch) and other collaborators.


License
=====

OpenMalaria is distributed under the terms of the
[GPL v2](http://opensource.org/licenses/GPL-2.0) (also see [COPYING](COPYING)),
or, (at your option) any later version.
