#!/bin/sh -e

if [ ! -d "vcpkg" ]; then
    git clone https://github.com/microsoft/vcpkg.git
fi

cd vcpkg && ./bootstrap-vcpkg.sh
cd ..
cp pyproject-vcpkg.toml pyproject.toml