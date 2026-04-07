import argparse
from pathlib import Path
from typing import Iterable
from sys import stdout

import numpy as np
from scipy import linalg
from loguru import logger as log

from . import _gdca
from . import _load_data

NAME = 'GaussDCApy'
VERSION = '1.0.0'
log.remove()
fmt = ('<green>{time:MM-DD HH:mm:ss}</green> | '
       '<level>{level: <8}</level> | '
       '<level>{message}</level>')
log.add(stdout, colorize=True, format=fmt, level='INFO', filter=NAME,
        backtrace=True, enqueue=True)
log.info(f'Starting {NAME.upper()} v{VERSION}')


def parse_fasta(fasta_file: Path) -> Iterable[tuple[str, str]]:
    """
    Import from treebarcode.
    Upper letter
    Return tuple(title: str, sequence: str)
    """
    with open(fasta_file, 'r') as f:
        record = []
        title = ''
        for line_raw in f:
            line = line_raw.strip()
            if line.startswith('>'):
                if len(record) != 0:
                    join_str = ''.join(record)
                    if len(join_str.strip()) != 0:
                        yield title, join_str.upper()
                    else:
                        log.warning(f'Found empty sequence {title}')
                    record.clear()
                title = line[1:]
            else:
                record.append(line)
        if len(record) != 0:
            join_str = ''.join(record)
            if len(join_str.strip()) != 0:
                yield title, join_str


def write_fasta(records: Iterable[tuple[str, str]], filename: Path) -> Path:
    # sequence is in one line, no truncated
    if filename.exists():
        log.warning(f'Overwriting existing file {filename}')
    with open(filename, 'w') as f:
        for title, seq in records:
            f.write(f'>{title}\n')
            f.write(f'{seq}\n')
    return filename


def array_to_fasta(name_array: np.ndarray, seq_array: np.ndarray,
                   output: Path) -> Path:
    # convert dtype='S1' to strings
    seq_list = list()
    for row in seq_array:
        seq_list.append(b''.join(row).decode('ascii'))
    return write_fasta(zip(name_array, seq_list), output)


def aln_to_array(records: Iterable[tuple[str, str]]) -> tuple[
    np.ndarray, np.ndarray]:
    records_ = list(records)
    name_array = np.array([i[0] for i in records_], dtype=np.str_)
    # S1 for 1 byte character, save more memory than 'U1' but
    # encode/decode is required
    seq_array = np.array(
        [np.fromiter(i[1], dtype=np.dtype('S1')) for i in records_])
    return name_array, seq_array


def _compute_FN(mJ, n_cols: int, alphabet_size: int):
    # int8 * int8 may overflow, which cause negative dimension in FN_all
    s = int(alphabet_size - 1)

    FN = np.zeros((n_cols, n_cols), dtype=np.float64)
    FN_all = np.zeros((n_cols, n_cols, s * s), dtype=np.float64)

    fs = s
    fs2 = s * s

    for i in range(n_cols - 1):
        _row = i * s
        for j in range(i + 1, n_cols):
            _col = j * s

            patch = mJ[_row: _row + s, _col: _col + s]
            total = patch.sum() / fs2
            rows = patch.sum(axis=1) / fs
            columns = patch.sum(axis=0) / fs

            fn_pre = patch - rows[:, None] - columns[None, :] + total
            fn = (fn_pre * fn_pre).sum()

            FN[i, j] = fn
            FN[j, i] = fn
            FN_all[i, j, :] = FN_all[j, i, :] = patch.flatten()

    FN = np.sqrt(FN)
    return FN, _gdca.apc_correction(FN), FN_all


def compute_ranking(scores, min_separation=5):
    # from MMichel/GaussDCA
    N = scores.shape[0]
    score_list = list()
    for i in range(N-min_separation):
        for j in range(i+min_separation, N):
            score_list.append((i+1, j+1, scores[i,j]))
    score_list.sort(key=lambda x:x[2], reverse=True)
    return score_list


def _compute_gdca_scores(alignment, alignment_T, min_separation=5):
    alphabet_size = alignment.max()

    n_cols = alignment_T.shape[1]
    depth = alignment_T.shape[0]

    covar, meff = _gdca.prepare_covariance(alignment, alignment_T)

    cho = linalg.cho_factor(covar, check_finite=False)
    mJ = linalg.cho_solve(cho, np.eye(covar.shape[0]), check_finite=False, overwrite_b=True)

    FN, FN_corr, FN_all = _compute_FN(mJ, n_cols, alphabet_size)
    covar_FN, covar_FN_corrected, covar_FN_all = _compute_FN(covar, n_cols, alphabet_size)
    results = dict(gdca=FN, gdca_corr=FN_corr, gdca_expanded=FN_all, eff_seq=meff, seq=depth,
                   covar_FN=covar_FN, covar_FN_corr=covar_FN_corrected, covar_expanded=covar_FN_all)
    score_list = compute_ranking(FN_corr, min_separation)
    return results, score_list


def run(align_file: str):
    output_file = Path(align_file).with_suffix('.txt')
    align = _load_data.load_fasta(align_file)
    results, score_list =  _compute_gdca_scores(np.ascontiguousarray(align), np.ascontiguousarray(align.T))
    with open(output_file, 'w') as f:
        for i, j, score in score_list:
            f.write(f'{i},{j},{score:.8f}\n')
    return score_list, output_file


def compute_weights(path, theta=None):
    ali = _load_data.load_fasta(path)
    if theta is None:
        theta = -1.

    return _gdca.compute_weights(np.ascontiguousarray(ali), np.ascontiguousarray(ali.T), theta)


def parse_args():
    arg = argparse.ArgumentParser()
    arg.add_argument('input', help='aligned protein sequence fasta file')
    arg.add_argument('-threshold', '-threshold', dest='threshold', default=0.75,
                     help='coevolution score threshold')
    return arg.parse_args()


def gaussdca_main():
    # Calculate DCA score and split the alignment
    arg = parse_args()
    fasta = Path(arg.input).resolve()
    assert fasta.exists()
    try:
        name, old_seq = aln_to_array(parse_fasta(fasta))
    except Exception:
        log.error(f'Bad alignment input {fasta}')
        raise SystemExit(-1)
    log.info(f'Input file: {fasta}')
    log.info(f'Coevolution score threshold: {arg.threshold}')
    output_file = fasta.with_suffix('.txt')
    co = fasta.with_suffix('.co.aln')
    non_co = fasta.with_suffix('.non_co.aln')
    co_half = fasta.with_suffix('.co_half.aln')
    non_co_half = fasta.with_suffix('.non_co_half.aln')
    invariant = fasta.with_suffix('.invariant.aln')
    all_except_half_co = fasta.with_suffix('.all_except_half_co.aln')

    unique_counts = np.array(
        [len(np.unique(old_seq[:, i])) for i in range(old_seq.shape[1])])
    invariant_index = np.where(unique_counts == 1)[0]
    invariant_site = old_seq[:, invariant_index]
    # mutant_index = np.where(unique_counts > 1)[0]

    score_list, score_txt = run(str(fasta))
    np_array = np.array(score_list)
    # fnr = call_julia(str(fasta))
    # # convert from numpy.void
    # np_array = np.array([list(i) for i in fnr.to_numpy()])
    # np.savetxt(result, np_array, fmt=['%d', '%d', '%.18e'])

    # start with 1->0
    np_array[:, :2] -= 1
    # strong coupling
    np_array2 = np_array[np_array[:, 2] > arg.threshold]
    all_index = set(np.arange(0, old_seq.shape[1]))
    co_index = set(np_array2[:, :2].flatten().astype(int))
    co_half_index = set(np_array2[:, 0].flatten().astype(int))
    non_co_index = all_index - co_index - set(invariant_index)
    non_co_half_index = all_index - co_half_index - set(invariant_index)
    all_except_half_index = all_index - co_half_index
    # output
    co_site = old_seq[:, list(co_index)]
    non_co_site = old_seq[:, list(non_co_index)]
    co_half_site = old_seq[:, list(co_half_index)]
    non_co_half_site = old_seq[:, list(non_co_half_index)]
    all_except_half_site = old_seq[:, list(all_except_half_index)]
    array_to_fasta(name, co_site, co)
    array_to_fasta(name, non_co_site, non_co)
    array_to_fasta(name, invariant_site, invariant)
    array_to_fasta(name, co_half_site, co_half)
    array_to_fasta(name, non_co_half_site, non_co_half)
    array_to_fasta(name, all_except_half_site, all_except_half_co)
    # print(result, co, non_co, co_half, non_co_half)
    log.info(f'{old_seq.shape[1]} columns')
    log.info(f'{np_array.shape[0]}, pairs')
    log.info(f'{np_array2.shape[0]}, pairs big score')
    log.info(f'{invariant_index.shape[0]} invariant sites')
    log.info(f'{len(co_index)} coevolution sites')
    log.info(f'{len(non_co_index)} non-coevolved sites')
    log.info(f'Output score file: {output_file}')
    log.info(f'Output alignments: {fasta.stem}.*.aln')
    log.info('Done!')
    return


if __name__ == '__main__':
    gaussdca_main()
