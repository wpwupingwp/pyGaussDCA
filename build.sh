python3 -m venv .venv
source .venv/bin/activate
pip install pythran scipy numpy wheel
pip install "setuptools<82"
python3 setup.py bdist_wheel
