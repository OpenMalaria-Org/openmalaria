#/bin/bash

# This file is part of OpenMalaria.
# 
# Copyright (C) 2005-2025 Swiss Tropical and Public Health Institute
# Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
# Copyright (C) 2020-2025 University of Basel
# Copyright (C) 2025 The Kids Research Institute Australia
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

# Assuming you have installed:
# gsl, git, cmake, xsd, xerces-c

# stop on error
set -e

# Clone to
OMGIT=openmalaria

# Clone from
RELEASE=openMalaria
BRANCH=main
SWITCHBRANCH=0

# options
CLEAN=0         # Clean build/
DEBUG=0
TESTS=OFF       # Don't generate tests
RUNTESTS=0      # Don't run tests
JOBS=4          # 4 threads
CREATERELEASE=0 # Create release ARTIFACT.zip:tar.gz
ARTIFACT=       # default filename is RELEASE-VERSION

# shell
CYGWIN=0

# For Windows - Cygwin (MobaXTerm)
export PATH=/usr/bin:$PATH

printHelp () {
    echo "HELP: make sure dependencies are installed (adapt to your distribution)"
    echo "======================================================================="
    echo "Ubuntu/Debian: sudo apt-get install build-essential git cmake libgsl-dev libxerces-c-dev xsdcxx"
    echo "--------------"
    echo "Mac OS: brew install git coreutils cmake gcc gsl xerces-c xsd"
    echo "--------------"
    echo "Windows - Cygwin (MobaXTerm): apt-get install p7zip gcc-g++ git cmake make zlib-devel libgsl-devel xsd libxerces-c-devel"

    echo ""
    echo "Options:"
    echo "  -b, --branch=<name>"    "specify the branch (master)"
    echo "  -c, --clean"            "clean build folder (false)"
    echo "  -d, --debug"            "build in debug mode (false)"
    echo "  -t, --tests"            "run the tests (false)"
    echo "  -r, --release"          "generate the release artifcat (false)"
    echo "  -a, --artifact=<name>"  "specify the artifcat name (openMalaria-VERSION)"
    echo "  -j, --jobs=<X>"         "specify the number of jobs (4)"
    echo "  -h, --help"             "print this message"
}

parseArguments () {
    for i in "$@"; do
        case $i in
            -b=*|--branch=*)    BRANCH="${i#*=}"; SWITCHBRANCH=1 && shift ;;
            -c|--clean)         CLEAN=1 && shift ;;
            -d|--debug)         DEBUG=1 && shift ;;
            -t|--tests)         TESTS=ON && shift ;;
            -r|--release)       CREATERELEASE=1 && shift ;;
            -a=*|--artifact=*)  ARTIFACT=${i#*=} && shift ;;
            -j=*|--jobs=*)      JOBS="${i#*=}" && shift ;;
            -h|--help)          printHelp && shift && exit ;;
            *) ;;
        esac
    done
}

isWindows () {
    unameOut="$(uname -s)"
    case "${unameOut}" in
        CYGWIN*)    CYGWIN=1;;
        MINGW*)     CYGWIN=1;;
        *)          CYGWIN=0;;
    esac
}

clone () {
    echo "Cloning..."
    # Already in openmalaria?
    if [ -d .git ] || [ git rev-parse --is-inside-work-tree ]; then
        if [ $(basename -s .git `git remote get-url origin`) = "openmalaria" ]; then
            echo "Already in openmalaria repo, not cloning."
        else
            echo "Error: this git repository is not openmalaria."
        fi
    else
        # Clone git repo
        if [ ! -d "$OMGIT" ] ; then
            git clone --branch $BRANCH https://github.com/OpenMalaria-Org/openmalaria.git $OMGIT
        else
            echo "Folder $OMGIT already exist, not cloning."
        fi
        cd $OMGIT
    fi
    if [ $SWITCHBRANCH -eq 1 ]; then
        echo "Switching branch to $BRANCH"
        git checkout $BRANCH && git pull
    fi
    git branch
}

build () {
    echo "Building..."
    mkdir -p build

    if [ $CLEAN -eq 1 ]; then
        cd build && rm -rf * && cd ..
    fi

    # Compile OpenMalaria
    cd build
    which cmake
    if [ $DEBUG -eq 1 ]; then
        cmake -DOM_BOXTEST_ENABLE=$TESTS -DOM_CXXTEST_ENABLE=$TESTS .. && make -j$JOBS
    else
        cmake -DCMAKE_BUILD_TYPE=Release -DOM_BOXTEST_ENABLE=$TESTS -DOM_CXXTEST_ENABLE=$TESTS .. && make -j$JOBS
    fi
    cd ..
}

runtests () {
    echo "Testing..."
    if [ $TESTS = "ON" ]; then
        cd build && ctest --output-on-failure -j$JOBS && cd ..
    fi
}

package () {
    echo "Packaging..."
    # Get version number
    VERSION=$(cat version.txt | cut -d'-' -f2)
    MAJOR=$(cat version.txt | cut -d'-' -f2 | cut -d'.' -f1)

    if [ -z "${ARTIFACT}" ]; then
        ARTIFACT=$RELEASE-$VERSION
        echo $ARTIFACT
    fi

    # Prepare release folder by copying the necessary files
    mkdir -p $ARTIFACT
    cp build/openMalaria $ARTIFACT/
    cp -r util/example/* $ARTIFACT/
    cp schema/scenario_$MAJOR.xsd $ARTIFACT/
    cp test/densities.csv $ARTIFACT/
    cp test/autoRegressionParameters.csv $ARTIFACT/
    cp COPYING $ARTIFACT/

    # if Cygwin, copy dll files
    if [ $CYGWIN -eq 1 ]; then
        cd $ARTIFACT/
        rm -f dlls
        for i in $(ldd openMalaria); do
            echo $i | grep "/usr" >> dlls || true
        done
        for i in $(cat dlls); do
            cp -f $i .
            echo "cp $i ."
        done
        rm -f dlls
        cd ..
    fi

    # Compress
    if [ $CYGWIN -eq 1 ]; then
        7z a $ARTIFACT.zip "$ARTIFACT/"
    else
        tar -czvf $ARTIFACT.tar.gz $ARTIFACT/
    fi
}

parseArguments $@ # CLEAN=0:1, JOBS=1:N, TESTS=ON:OFF, print help
isWindows # set CYGWIN=0:1
clone # clone the repo
build # compile

if [ $TESTS = "ON" ]; then
    runtests # test
fi

if [ $CREATERELEASE -eq 1 ]; then
    package # create $ARTIFACT.zip or $ARTIFACT.tar.gz
    echo "RELEASE: $(ls $ARTIFACT.*)"
fi

# Don't delete the temp folders yet, let the user decide
# rm -rf $OMGIT $ARTIFACT/