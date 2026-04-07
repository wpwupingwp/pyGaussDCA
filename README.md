# pyGaussDCA2
Reimplementation of [pyGaussDCA](https://github.com/ElofssonLab/pyGaussDCA)

Used output function from https://github.com/MMichel/GaussDCA.

## New features
1. add pyproject.toml
2. add txt output (like GaussDCA.jl)
3. robust fasta parser
4. add wheel files for Python 3.10-3.14

## Usage
```bash
pip install pyGaussDCA2
python3 -m GaussDCA input.fasta
```

## Compatibility
Currently only support Linux due to pythran compile failed in Windows


