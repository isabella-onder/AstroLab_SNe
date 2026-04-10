import numpy as np
from pathlib import Path

COLOUR_PAIRS = [
    ("Bs", "Vs"),
    #("Vs", "Rs"),
    ("Rs", "Is"),
    ("Bs", "Is"),
    ("Vs", "Is"),
]

MAX_ERR = 0.7

def parse_file(filepath):
    meta, data, current_filter = {}, {}, None
    mjd_buf, mag_buf, err_buf = [], [], []

    def _flush():
        if current_filter is not None and mjd_buf:
            data[current_filter] = {
                "mjd":     np.array(mjd_buf,  dtype=float),
                "mag":     np.array(mag_buf,  dtype=float),
                "mag_err": np.array(err_buf,  dtype=float),
            }

    with open(filepath) as fh:
        lines = [ln.strip() for ln in fh if ln.strip()]

    header = lines[0].split()
    meta = {"name": header[0], "redshift": float(header[1]),
            "ra": float(header[2]), "dec": float(header[3])}

    for line in lines[1:]:
        parts = line.split()
        if parts[0].lower() == "filter":
            _flush()
            current_filter = parts[1]
            mjd_buf, mag_buf, err_buf = [], [], []
        elif len(parts) == 3:
            try:
                mjd_buf.append(float(parts[0]))
                mag_buf.append(float(parts[1]))
                err_buf.append(float(parts[2]))
            except ValueError:
                pass
    _flush()
    return meta, {k: v for k, v in data.items() if v is not None}


def match_epochs(d1, d2, tol=0.3):
    mjd_out, col_out, err_out = [], [], []
    for i, t1 in enumerate(d1["mjd"]):
        dt = np.abs(d2["mjd"] - t1)
        match = np.where(dt <= tol)[0]
        if not match.size:
            continue
        j = match[np.argmin(dt[match])]
        mjd_out.append(0.5 * (t1 + d2["mjd"][j]))
        col_out.append(d1["mag"][i] - d2["mag"][j])
        err_out.append(np.sqrt(d1["mag_err"][i]**2 + d2["mag_err"][j]**2))
    return np.array(mjd_out), np.array(col_out), np.array(err_out)


def get_colour_curves(sn_name, filepath=None, colour_pairs=COLOUR_PAIRS):
    """
    Parameters
    ----------
    sn_name : str
        SN name used to locate the file if filepath not given.
    filepath : str or Path, optional
        Explicit path to the photometry file.

    Returns
    -------
    list of [label, mjd, colour, colour_err]
        One entry per available colour pair, where:
          label      : str   e.g. 'Bs-Vs'
          mjd        : array of epoch midpoints
          colour     : array of mag differences (filter1 - filter2)
          colour_err : array of propagated uncertainties
    """
    fp = Path(filepath) if filepath else Path(f"{sn_name}.txt")
    _, data = parse_file(fp)

    def quality_cut(d):
        mask = d["mag_err"] <= MAX_ERR
        return {k: v[mask] for k, v in d.items()}

    results = []
    for f1, f2 in colour_pairs:
        if f1 not in data or f2 not in data:
            continue
        mjd, col, err = match_epochs(quality_cut(data[f1]), quality_cut(data[f2]))
        if mjd.size:
            results.append([f"{f1}-{f2}", mjd, col, err])

    return results