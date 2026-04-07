#!/bin/bash
version_list=("3.10" "3.11" "3.12" "3.13" "3.14")
for version in ${version_list[@]}
do
    rm -rf build
    mkdir build
    uv python pin $version
    uv sync
    rm -rf src/gaussdcapy/*.cpp
    rm -rf src/*info
    rm -rf build/*
    rm -rf dist/gaussdcapy
    rm -rf dist/*info
    #uv run python3 setup.py build
    #uv run python3 setup.py build_ext
    uv run python3 setup.py bdist_wheel
done
