Open Malaria
============

This is a pre-built version of OpenMalaria. For documentation, see:
https://github.com/OpenMalaria-Org/openmalaria/wiki

### What's included ###

```
openMalaria / openMalaria.exe       — the executable
lib*.so                             — unix library files
*.dll                               — windows library files
*.csv                               — resource files
scenario_*.xsd                      — schema file

example_scenario.xml                — an example
run-example-scenario.*              — script to run the included example
output.txt, ctsout.txt              — output from this example

version.txt                         — the OpenMalaria version
appveyor*, travis*                  — build info
, dependencies.txt
```

Library files (.so, .dll), the schema (.xsd) and resource files (.csv) should
be in the same directory as the executable when running OpenMalaria. Not all
may be required, depending on your system libraries and resources needed by the
simulation. Alternatively resource files may be retrieved from a different path
via the --resource-path option.

Linux builds are currently missing required libraries. If your system has
different library versions, it is probably easiest to clone the repository and
build from source.

Windows builds currently include all schema versions in the `schema/` directory.
Copy the correct version for this OpenMalaria build to the current directory.
(This must correspond to the build version, not the scenario version.)

### Scenario version ###

Each scenario starts similar to the following:
```xml
<?xml version='1.0' encoding='UTF-8'?>
<om:scenario xmlns:om="http://openmalaria.org/schema/scenario_41" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" name="Example Scenario" schemaVersion="41" xsi:schemaLocation="http://openmalaria.org/schema/scenario_41 scenario_41.xsd">
```

In case the version number does not match the version of this OpenMalaria build,
it must be updated to match (the major part of the version number only). Here:

-   `schemaVersion="41"` — the schema version is 41
-   `http://openmalaria.org/schema/scenario_41` — this is the XML namespace
-   `scenario_41.xsd` — the name of the schema file

For the most part, this is the only change needed when updating an XML to use a
newer version of OpenMalaria, though this is not always the case.


## License ##

OpenMalaria is distributed under the terms of the
[GPL v2](http://opensource.org/licenses/GPL-2.0) (also see [COPYING](COPYING)),
or, (at your option) any later version.
