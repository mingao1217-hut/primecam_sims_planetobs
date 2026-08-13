#!/usr/bin/env python3
"""
Per-detector minimum distance from each true detector position to the
actual simulated Jupiter trajectory, for one HDF5 dataset -- identical
computation to diagnose_detector_proximity.py, but additionally broken
down per wafer_slot (w1/w2/w3), since a full 3-array module truth
focalplane may have very different coverage per wafer depending on scan
geometry and hardware tube offset.

Kept as a separate file from diagnose_detector_proximity.py (which has no
wafer_slot awareness at all) rather than modifying it, so the existing
single-array (w2-only) diagnostic and all its prior outputs stay exactly
reproducible.

For every detector i in the full truth focal plane:
    dmin_i = min_t hypot(x_traj(t) - x_det_i, y_traj(t) - y_det_i)

No detector timestream data is used -- purely boresight pointing (times +
boresight_azel) against the truth focal plane, so this works even when
only a handful of detectors were actually simulated (as in the stage-7
10-det trajectory-generation run): dmin is computed against every
detector listed in --truth-focalplane, independent of how many were
simulated.

Example
-------
python diagnose_detector_proximity_wafer.py \
  --input-dir ./ccat_datacenter_mock/mockdata/planet_ATMdata_d10_full_module_i6_d10/sim_PCAM280_h5_Jupiter_d10/ \
  --truth-focalplane input_files/fp_files/fp_f280_dettable_I6.h5 \
  --truth-table dettable_stack \
  --out-prefix proximity_full_module_i6_d10
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import toast
from astropy.table import QTable

from fit_all_detectors_beam import (
    azel_to_toast_vector,
    jupiter_azel,
    pointing_vector_from_quat,
    quat_to_matrix_xyzw,
    tangent_xy_arcmin,
    text_value,
)


def load_truth_with_wafer(path, table_path):
    table = QTable.read(path, path=table_path)
    for col in ("name", "quat", "wafer_slot"):
        if col not in table.colnames:
            raise KeyError(f"{path}:{table_path} must contain a '{col}' column.")

    names = [text_value(n).strip() for n in table["name"]]
    wafer_slot = np.array([text_value(w).strip() for w in table["wafer_slot"]])
    quats = np.asarray(table["quat"], dtype=float)
    vec = pointing_vector_from_quat(quats)
    x, y = tangent_xy_arcmin(vec)

    return names, wafer_slot, x, y


def print_stats(label, dmin, fwhm_arcmin):
    half_fwhm = 0.5 * fwhm_arcmin
    one_fwhm = fwhm_arcmin
    print(f"{label} (N={len(dmin)}):")
    print(f"    median: {np.median(dmin):.3f}'")
    print(f"    p90:    {np.percentile(dmin, 90):.3f}'")
    print(f"    p95:    {np.percentile(dmin, 95):.3f}'")
    print(f"    max:    {np.max(dmin):.3f}'")
    print(f"    fraction with dmin <= 0.5 FWHM ({half_fwhm:.4f}'): "
          f"{np.mean(dmin <= half_fwhm):.4f}")
    print(f"    fraction with dmin <= 1.0 FWHM ({one_fwhm:.4f}'):  "
          f"{np.mean(dmin <= one_fwhm):.4f}")
    print()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--truth-focalplane", type=Path, required=True)
    parser.add_argument("--truth-table", default="dettable_stack")
    parser.add_argument("--fwhm-arcmin", type=float, default=0.7833)
    parser.add_argument("--det-batch-size", type=int, default=200,
                         help="Detectors processed per chunk, to bound memory "
                              "(full N_det x N_samples distance matrix is avoided)")
    parser.add_argument("--color-max-arcmin", type=float, default=5.0,
                         help="Fixed colorbar max for the focal-plane plots, so "
                              "different passes are directly comparable")
    parser.add_argument("--out-prefix", default="detector_proximity_wafer")
    args = parser.parse_args()

    world, procs, rank = toast.get_world()
    if procs != 1:
        raise RuntimeError("Run with plain `python`, not mpirun.")

    data = toast.Data(comm=toast.Comm(world=world, groupsize=1))

    loader = toast.ops.LoadHDF5(
        volume=str(args.input_dir),
        pattern=r".*\.h5$",
        meta=["name", "uid", "telescope", "session"],
        shared=["times", "boresight_azel"],
        detdata=[],
        sort_by_size=False,
        process_rows=1,
        force_serial=True,
    )
    loader.apply(data)

    if len(data.obs) == 0:
        raise RuntimeError(f"No observations found in {args.input_dir}")
    if len(data.obs) > 1:
        raise RuntimeError(
            f"{len(data.obs)} observations found -- this script assumes a single, "
            "single-scan dataset."
        )
    obs = data.obs[0]

    det_names, wafer_slot, x_det, y_det = load_truth_with_wafer(
        args.truth_focalplane, args.truth_table
    )
    n_det = len(det_names)

    times = np.asarray(obs.shared["times"].data, dtype=float)
    bore_quat = np.asarray(obs.shared["boresight_azel"][:], dtype=float)

    jup_az, jup_el = jupiter_azel(times, obs.telescope.site)
    jup_world = azel_to_toast_vector(jup_az, jup_el)
    R_bore = quat_to_matrix_xyzw(bore_quat)
    jup_local = np.einsum("...ji,...j->...i", R_bore, jup_world)
    x_traj, y_traj = tangent_xy_arcmin(jup_local)

    print(f"{args.input_dir}")
    print(f"{obs.name}: N_samples={len(times)}, N_det={n_det} "
          f"(truth focalplane: {args.truth_focalplane}:{args.truth_table})\n")

    dmin = np.full(n_det, np.inf)
    for start in range(0, n_det, args.det_batch_size):
        stop = min(start + args.det_batch_size, n_det)
        dx = x_det[start:stop, None] - x_traj[None, :]
        dy = y_det[start:stop, None] - y_traj[None, :]
        dist = np.hypot(dx, dy)
        dmin[start:stop] = dist.min(axis=1)

    print_stats("Overall", dmin, args.fwhm_arcmin)

    wafers = sorted(np.unique(wafer_slot))
    for w in wafers:
        mask = wafer_slot == w
        print_stats(f"Wafer {w}", dmin[mask], args.fwhm_arcmin)

    one_fwhm = args.fwhm_arcmin

    # Overall plot (same style as diagnose_detector_proximity.py)
    fig, ax = plt.subplots(figsize=(8, 10))
    sc = ax.scatter(
        x_det, y_det, c=dmin, cmap="viridis_r", s=8,
        vmin=0.0, vmax=args.color_max_arcmin,
    )
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04, extend="max")
    cbar.set_label("Minimum distance to Jupiter trajectory [arcmin]")
    ax.set_xlabel("dxel [arcmin]")
    ax.set_ylabel("del [arcmin]")
    ax.set_title(
        f"{args.out_prefix}: per-detector min. distance to Jupiter trajectory\n"
        f"(median={np.median(dmin):.2f}', frac<=1 FWHM={np.mean(dmin <= one_fwhm):.1%})"
    )
    ax.grid(alpha=0.3)
    ax.set_aspect("equal")
    fig.tight_layout()
    out_png = Path(f"{args.out_prefix}.png")
    fig.savefig(out_png, dpi=150)
    print(f"Wrote: {out_png}")

    # Per-wafer breakdown plot: one panel per wafer, same color scale
    fig2, axes = plt.subplots(1, len(wafers), figsize=(6 * len(wafers), 7), squeeze=False)
    axes = axes[0]
    for ax, w in zip(axes, wafers):
        mask = wafer_slot == w
        sc = ax.scatter(
            x_det[mask], y_det[mask], c=dmin[mask], cmap="viridis_r", s=10,
            vmin=0.0, vmax=args.color_max_arcmin,
        )
        ax.set_xlabel("dxel [arcmin]")
        ax.set_ylabel("del [arcmin]")
        ax.set_aspect("equal")
        ax.grid(alpha=0.3)
        ax.set_title(
            f"{w} (N={mask.sum()})\n"
            f"median={np.median(dmin[mask]):.2f}', "
            f"frac<=1 FWHM={np.mean(dmin[mask] <= one_fwhm):.1%}"
        )
    cbar = fig2.colorbar(sc, ax=axes, fraction=0.02, pad=0.02, extend="max")
    cbar.set_label("Minimum distance to Jupiter trajectory [arcmin]")
    fig2.suptitle(f"{args.out_prefix}: per-wafer breakdown")
    out_png_wafer = Path(f"{args.out_prefix}_by_wafer.png")
    fig2.savefig(out_png_wafer, dpi=150)
    print(f"Wrote: {out_png_wafer}")


if __name__ == "__main__":
    main()
