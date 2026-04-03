from pathlib import Path
from sys import argv

import numpy as np
from scipy import linalg
from operator import itemgetter

from . import _gdca
from . import _load_data


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
    N = scores.shape[0]
    score_list = list()
    for i in range(N-min_separation):
        for j in range(i+min_separation, N):
            score_list.append((i+1, j+1, scores[i,j]))
    score_list.sort(key=itemgetter(2), reverse=True)
    return score_list


def _compute_gdca_scores(alignment, alignment_T, verbose, min_separation=5):
    alphabet_size = alignment.max()

    n_cols = alignment_T.shape[1]
    depth = alignment_T.shape[0]

    if verbose:
        print('Prepare inputs')
    covar, meff = _gdca.prepare_covariance(alignment, alignment_T)

    if verbose:
        print('Invert matrix')
    cho = linalg.cho_factor(covar, check_finite=False)
    mJ = linalg.cho_solve(cho, np.eye(covar.shape[0]), check_finite=False, overwrite_b=True)

    if verbose:
        print('Compute Frobenius Norm')
    FN, FN_corr, FN_all = _compute_FN(mJ, n_cols, alphabet_size)
    covar_FN, covar_FN_corrected, covar_FN_all = _compute_FN(covar, n_cols, alphabet_size)
    results = dict(gdca=FN, gdca_corr=FN_corr, gdca_expanded=FN_all, eff_seq=meff, seq=depth,
                   covar_FN=covar_FN, covar_FN_corr=covar_FN_corrected, covar_expanded=covar_FN_all)
    score_list = compute_ranking(FN_corr, min_separation)
    return results, score_list


def run(align_file: str, verbose=False):
    if verbose:
        print('Loading data')
    output_file = Path(align_file).with_suffix('.txt')
    align = _load_data.load_a3m(align_file)
    results, score_list =  _compute_gdca_scores(np.ascontiguousarray(align), np.ascontiguousarray(align.T), verbose)
    with open(output_file, 'w') as f:
        for i, j, score in score_list:
            f.write(f'{i},{j},{score}\n')
    print(f'Output file: {output_file}')
    return output_file


def compute_weights(path, theta=None):
    ali = _load_data.load_a3m(path)
    if theta is None:
        theta = -1.

    return _gdca.compute_weights(np.ascontiguousarray(ali), np.ascontiguousarray(ali.T), theta)


def main():
    # aligned a3m or fasta (one line)
    align_file = Path(argv[1]).resolve()
    assert align_file.exists()
    run(str(align_file), verbose=False)
    return


if __name__ == '__main__':
    main()
