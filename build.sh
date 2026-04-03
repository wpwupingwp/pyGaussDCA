#!/bin/bash
uv run python3 setup.py build
uv run python3 setup.py build_ext
uv run python3 setup.py bdist_wheel
