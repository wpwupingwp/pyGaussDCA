# GaussDCApy
Reimplementation of [pyGaussDCA](https://github.com/ElofssonLab/pyGaussDCA)

Used output function from https://github.com/MMichel/GaussDCA.

## New features
1. add pyproject.toml
2. add txt output (like GaussDCA.jl)
3. use robust fasta parser instead of a3m parser
4. add wheel files for Python 3.10-3.14
5. split alignment according to coevolution score

## Usage
```bash
pip install gaussdcapy
python3 -m gaussdcapy input.fasta
```
If there is no valid wheel files, see `build.sh` for building the package from
source.

## Compatibility
Currently only support Linux due to pythran compile failed in Windows.
