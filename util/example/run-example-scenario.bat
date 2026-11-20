REM This file is part of OpenMalaria.
REM 
REM Copyright (C) 2005-2026 Swiss Tropical and Public Health Institute
REM Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
REM Copyright (C) 2020-2026 University of Basel
REM Copyright (C) 2025-2026 The Kids Research Institute Australia
REM
REM OpenMalaria is free software; you can redistribute it and/or modify
REM it under the terms of the GNU General Public License as published by
REM the Free Software Foundation; either version 2 of the License, or (at
REM your option) any later version.
REM 
REM This program is distributed in the hope that it will be useful, but
REM WITHOUT ANY WARRANTY; without even the implied warranty of
REM MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
REM General Public License for more details.
REM 
REM You should have received a copy of the GNU General Public License
REM along with this program; if not, write to the Free Software
REM Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

del output.txt
del ctsout.txt
openMalaria.exe --scenario example_scenario.xml
pause
