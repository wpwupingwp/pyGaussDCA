#!/bin/bash
rm -rf src/gaussdca/*.cpp
rm -rf build/*
rm -rf dist/gaussdca
rm -rf dist/pygaussdca-0.1.0.dist-info
rm -rf dist/pygaussdca-0.1.0-cp312-cp312-linux_x86_64.whl
#uv run python3 setup.py build
#uv run python3 setup.py build_ext
uv run python3 setup.py bdist_wheel
