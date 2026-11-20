Open Malaria
============

This is a pre-built version of OpenMalaria. For documentation, see:
https://github.com/OpenMalaria-Org/openmalaria/wiki

### What's included

```
openMalaria / openMalaria.exe       — the executable
*.dll                               — library files (Windows only)
*.csv                               — resource files
scenario_*.xsd                      — schema file
example_scenario.xml                — an example scenario
run-example-scenario.bat            — script to run the example on Windows
run-example-scenario.sh             — script to run the example on Linux and MacOS
COPYING                             - a copy of the GNU GPL v2 license text used by OpenMalaria
```

Library files (.dll), the schema (.xsd) and resource files (.csv) should
be in the same directory as the executable when running OpenMalaria. Not all
may be required, depending on your system libraries and resources needed by the
simulation. Alternatively resource files may be retrieved from a different path
via the --resource-path option.

Linux builds target specific Linux distribution and versions. If your system has
different library versions, it is probably easiest to clone the repository and
build from source, see: https://github.com/OpenMalaria-Org/openmalaria/wiki/BuildingOpenMalaria.

### Scenario version

Each scenario starts similar to the following:
```xml
<?xml version='1.0' encoding='UTF-8'?>
<om:scenario xmlns:om="http://openmalaria.org/schema/scenario_50" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" name="Example Scenario" schemaVersion="50" xsi:schemaLocation="http://openmalaria.org/schema/scenario_50 scenario_50.xsd">
```

In case the version number does not match the version of this OpenMalaria build,
it must be updated to match (the major part of the version number only). Here:

-   `schemaVersion="50"` — the schema version is 50
-   `http://openmalaria.org/schema/scenario_50` — this is the XML namespace
-   `scenario_50.xsd` — the name of the schema file

For the most part, this is the only change needed when updating an XML to use a
newer version of OpenMalaria, though this is not always the case, see:
https://github.com/OpenMalaria-Org/openmalaria/wiki/SchemaUpdateGuide
https://github.com/OpenMalaria-Org/openmalaria/wiki/Changelog


## License

OpenMalaria is distributed under the terms of the GNU General Public License version 2 (GPL v2).
A copy of the license is included in this folder, see COPYING.
Also see: http://opensource.org/licenses/GPL-2.0.
