#!/bin/bash

# This file is part of OpenMalaria.
# 
# Copyright (C) 2005-2026 Swiss Tropical and Public Health Institute
# Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
# Copyright (C) 2020-2026 University of Basel
# Copyright (C) 2025-2026 The Kids Research Institute Australia
#
# OpenMalaria is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

# Profile-guided optimisation (PGO) build of openMalaria with GCC or Clang.
# Builds an instrumented binary, runs a few test scenarios to collect a
# profile, then rebuilds using that profile. Output is identical to a normal
# Release build; only speed differs. Used by `build.sh --pgo`.
#
# usage: util/pgo-build.sh [build-dir] [extra cmake args]
# JOBS sets make parallelism (default: number of CPUs).

set -euo pipefail
SRC=$(cd "$(dirname "$0")/.." && pwd)
OUT=${1:-$SRC/build-pgo}; shift || true
mkdir -p "$OUT" && OUT=$(cd "$OUT" && pwd)
PROF=$OUT/pgo-profile
JOBS=${JOBS:-$(getconf _NPROCESSORS_ONLN)}
# Mix of within-host and transmission models; add scenarios close to your own workload.
TRAIN="Genotypes Penny VecFullTest NonHumanHosts MSAT_HRP2 Vivax"

if "${CXX:-c++}" --version | grep -qi clang; then
  GEN="-fprofile-generate=$PROF"
  USE="-fprofile-use=$PROF/merged.profdata -Wno-profile-instr-unprofiled -Wno-profile-instr-out-of-date"
  PROFDATA=$(command -v llvm-profdata || xcrun -f llvm-profdata)
else
  # The same build dir is reused for both phases: gcc keys profiles by object path.
  GEN="-fprofile-generate -fprofile-update=single -fprofile-dir=$PROF"
  USE="-fprofile-use -fprofile-partial-training -fprofile-dir=$PROF -Wno-missing-profile"
fi

rm -rf "$PROF" && mkdir -p "$PROF"

# 1. Instrumented build.
cmake -S "$SRC" -B "$OUT" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="$GEN" "$@"
make -C "$OUT" -j"$JOBS"

# 2. Training runs.
for s in $TRAIN; do
  d=$(mktemp -d)
  cp "$OUT"/schema/scenario_current.xsd "$SRC"/schema/*.xsd "$d"/
  (cd "$d" && "$OUT/openMalaria" --resource-path "$SRC/test" --scenario "$SRC/test/scenario$s.xml" >/dev/null)
  rm -rf "$d"
done
[ -n "${PROFDATA:-}" ] && "$PROFDATA" merge -output="$PROF/merged.profdata" "$PROF"/*.profraw

# 3. Optimised build using the profile. CMAKE_CXX_FLAGS stays in the cache, so
# later reconfigures of this build dir (e.g. build.sh enabling tests) keep PGO.
make -C "$OUT" clean
cmake -S "$SRC" -B "$OUT" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="$USE" "$@"
make -C "$OUT" -j"$JOBS"
echo "PGO build: $OUT/openMalaria"
