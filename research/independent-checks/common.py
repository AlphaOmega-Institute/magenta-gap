"""Shared input loader for the independent checks.

Inputs   research/object-colour/inputs/stockman_sharpe_2deg_390_830.csv: the 441 Stockman-Sharpe 2-degree LMS samples,
         390-830 nm at 1 nm, as exported for the continuum rim certificate (see the licence file beside it).
Returns  wavelengths (441,) and LMS (441, 3) as float64, exactly as stored.
"""
import csv
import pathlib

import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[2]
TABLE = REPO / "research" / "object-colour" / "inputs" / "stockman_sharpe_2deg_390_830.csv"


def load_table():
    rows = [r for r in csv.reader(TABLE.read_text().splitlines()) if r]
    data = np.array([[float(v) for v in r] for r in rows[1:]])
    assert data.shape == (441, 4) and data[0, 0] == 390.0 and data[-1, 0] == 830.0
    return data[:, 0].copy(), data[:, 1:].copy()
