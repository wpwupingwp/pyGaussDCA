#!/bin/bash
rm -rf src/gaussdca/*.cpp
rm -rf build/*
rm -rf dist/gaussdca
rm -rf dist/pygaussdca-*.dist-info
#uv run python3 setup.py build
#uv run python3 setup.py build_ext
uv run python3 setup.py bdist_wheel
