import csv
import re
from collections import Counter
from pathlib import Path
from statistics import median

import numpy as np
import wfdb
from wfdb import processing
from scipy.signal import butter, sosfiltfilt

WORK = Path.home() / "Downloads/POTASSIUM_ECG_RESEARCH"
ROOT = WORK / "waveforms_mimic_iv_ecg_v1"

QC = WORK / "potassium_waveform_QC_PRIVATE.csv"
PAIRS = WORK / "potassium_serial_pairs_PRIVATE.csv"

FEATURES = WORK / "potassium_geometry_features_v2_PRIVATE.csv"
RESULTS = WORK / "potassium_geometry_pair_results_v2_PRIVATE.csv"
MATRICES = WORK / "potassium_geometry_matrices_v2_PRIVATE.npz"

LEADS = ["I", "II", "V1", "V2", "V3", "V4", "V5", "V6"]

for p in (QC, PAIRS):
    if not p.exists():
        raise SystemExit(f"Missing input: {p}")

for p in (FEATURES, RESULTS, MATRICES):
    if p.exists():
        raise SystemExit(f"Output already exists: {p.name}")

def normalize(value):
    s = str(value).strip().replace("\\", "/")
    s = re.sub(r"\.(hea|dat)$", "", s)
    m = re.search(
        r"(files/p\d{4}/p\d{8}/s\d{8}/\d{8})$", s
    )
    if not m:
        raise ValueError("Unexpected waveform path")
    return m.group(1)

def repair_small_gaps(x):
    x = x.copy()

    for j in range(x.shape[1]):
        y = x[:, j]
        good = np.isfinite(y)

        if good.all():
            continue

        bad = ~good

        if bad.sum() > 25 or not good[0] or not good[-1]:
            raise ValueError("Excessive missing samples")

        transitions = np.diff(
            np.r_[False, bad, False].astype(int)
        )
        starts = np.flatnonzero(transitions == 1)
        ends = np.flatnonzero(transitions == -1)

        if max(ends - starts, default=0) > 5:
            raise ValueError("Long missing-signal gap")

        indices = np.arange(len(y))

        x[:, j] = np.interp(
            indices, indices[good], y[good]
        )

    return x

def geometry(matrix):
    centered = matrix - matrix.mean(
        axis=1, keepdims=True
    )

    u, s, _ = np.linalg.svd(
        centered, full_matrices=False
    )

    energy = s ** 2
    total = energy.sum()

    if total < 1e-12:
        raise ValueError("Insufficient matrix energy")

    return (
        float(energy[0] / total),
        float(energy[3:].sum() / total),
        u[:, 0]
    )

def axis_angle(a, b):
    alignment = abs(float(np.dot(a, b)))
    alignment = np.clip(alignment, 0, 1)
    return float(np.degrees(np.arccos(alignment)))

def extract(path):
    record = wfdb.rdrecord(str(ROOT / path))
    x = np.asarray(record.p_signal, dtype=float)

    if x.ndim != 2 or x.shape[0] != 5000:
        raise ValueError("Unexpected waveform dimensions")

    if abs(float(record.fs) - 500) > 0.01:
        raise ValueError("Unexpected sampling frequency")

    names = list(record.sig_name)

    indices = [names.index(lead) for lead in LEADS]
    x = repair_small_gaps(x[:, indices])

    fs = 500

    # Exploratory morphology filtering.
    sos = butter(
        3, [0.05, 40],
        btype="bandpass",
        fs=fs,
        output="sos"
    )

    y = sosfiltfilt(sos, x, axis=0)
    best = None

    for lead in ("II", "V2", "V4", "I", "V5", "V1"):
        channel = LEADS.index(lead)

        try:
            peaks = processing.xqrs_detect(
                sig=y[:, channel],
                fs=fs,
                verbose=False
            )
        except Exception:
            continue

        peaks = np.asarray(peaks, dtype=int)

        if not 4 <= len(peaks) <= 35:
            continue

        rr = np.diff(peaks) / fs

        if not 0.30 <= float(np.median(rr)) <= 2.0:
            continue

        qrs_beats = []
        post_beats = []
        durations = []

        for j, peak in enumerate(peaks[:-1]):
            if peak - 75 < 0 or peak + 50 > len(y):
                continue

            next_peak = peaks[j + 1]

            # QRS: fixed 200-ms window.
            qrs_start = peak - 50
            qrs_end = peak + 50

            # Adaptive post-QRS segment.
            post_start = peak + 70
            post_end = min(
                peak + 220,
                next_peak - 40
            )

            # Minimum 100 ms of post-QRS signal.
            if post_end - post_start < 50:
                continue

            baseline = np.median(
                y[peak-75:peak-55, :],
                axis=0
            )

            qrs = (
                y[qrs_start:qrs_end, :] - baseline
            ).T

            segment = (
                y[post_start:post_end, :] - baseline
            ).T

            # Normalize the variable-length segment
            # to 150 points for geometric comparison.
            original_time = np.linspace(
                0, 1, segment.shape[1]
            )
            target_time = np.linspace(0, 1, 150)

            post = np.stack([
                np.interp(
                    target_time, original_time, lead_signal
                )
                for lead_signal in segment
            ])

            qrs_beats.append(qrs)
            post_beats.append(post)

            durations.append(
                (post_end - post_start) / fs * 1000
            )

        if len(qrs_beats) < 2:
            continue

        qrs_beats = np.asarray(qrs_beats)
        post_beats = np.asarray(post_beats)

        reference = np.median(
            qrs_beats, axis=0
        ).ravel()

        reference = reference - reference.mean()
        ref_norm = np.linalg.norm(reference)

        correlations = []

        for beat in qrs_beats:
            b = beat.ravel()
            b = b - b.mean()

            denom = np.linalg.norm(b) * ref_norm

            corr = (
                float(np.dot(b, reference) / denom)
                if denom > 1e-10 else 0.0
            )

            correlations.append(corr)

        accepted = np.array(correlations) >= 0.75

        if accepted.sum() < 2:
            continue

        score = (
            int(accepted.sum()),
            float(np.median(
                np.asarray(correlations)[accepted]
            ))
        )

        if best is None or score > best[0]:
            best = (
                score,
                np.median(
                    qrs_beats[accepted], axis=0
                ),
                np.median(
                    post_beats[accepted], axis=0
                ),
                float(np.median(
                    np.asarray(durations)[accepted]
                )),
                float(np.median(rr) * 1000),
                lead
            )

    if best is None:
        raise ValueError("No consistent adaptive beats")

    score, qrs, post, post_ms, rr_ms, lead = best

    qpc1, qtail, qaxis = geometry(qrs)
    ppc1, ptail, paxis = geometry(post)

    features = {
        "qrs_rms": float(
            np.sqrt(np.mean(qrs ** 2))
        ),
        "post_rms": float(
            np.sqrt(np.mean(post ** 2))
        ),
        "qrs_pc1": qpc1,
        "post_pc1": ppc1,
        "qrs_tail3": qtail,
        "post_tail3": ptail,
        "qrs_post_angle": axis_angle(
            qaxis, paxis
        ),
        "post_trajectory_length": float(
            np.linalg.norm(
                np.diff(post, axis=1), axis=0
            ).sum()
        ),
        "median_rr_ms": rr_ms,
        "post_window_ms": post_ms,
        "beats_used": score[0],
        "beat_consistency": score[1],
        "detector_lead": lead
    }

    return (
        features,
        qrs.astype("float32"),
        post.astype("float32")
    )

with QC.open(encoding="utf-8-sig", newline="") as f:
    paths = sorted(
        normalize(row["waveform"])
        for row in csv.DictReader(f)
        if row["status"] == "OK"
    )

print("\n===== ADAPTIVE ECG GEOMETRY V2 =====")
print("Clean ECGs to process:", len(paths), flush=True)

features = {}
matrices = {}
failures = Counter()

for i, path in enumerate(paths, 1):
    try:
        feat, qrs, post = extract(path)
        features[path] = feat
        matrices[path] = (qrs, post)

    except Exception as exc:
        failures[
            type(exc).__name__ + ": " + str(exc)
        ] += 1

    if i % 100 == 0 or i == len(paths):
        print(
            f"Processed {i}/{len(paths)}",
            flush=True
        )

columns = [
    "waveform",
    "qrs_rms",
    "post_rms",
    "qrs_pc1",
    "post_pc1",
    "qrs_tail3",
    "post_tail3",
    "qrs_post_angle",
    "post_trajectory_length",
    "median_rr_ms",
    "post_window_ms",
    "beats_used",
    "beat_consistency",
    "detector_lead"
]

with FEATURES.open(
    "x", encoding="utf-8-sig", newline=""
) as f:
    writer = csv.DictWriter(
        f, fieldnames=columns
    )
    writer.writeheader()

    for path in sorted(features):
        writer.writerow({
            "waveform": path,
            **features[path]
        })

order = sorted(matrices)

if not order:
    raise SystemExit("No matrices extracted.")

np.savez_compressed(
    MATRICES,
    waveform=np.array(order),
    qrs=np.stack([
        matrices[p][0] for p in order
    ]),
    post_qrs=np.stack([
        matrices[p][1] for p in order
    ])
)

pair_results = []
counts = Counter()

with PAIRS.open(encoding="utf-8-sig", newline="") as f:
    for row in csv.DictReader(f):
        group = row["group"]
        counts[group] += 1

        a = normalize(row["abnormal_waveform"])
        n = normalize(row["normal_waveform"])

        if a not in matrices or n not in matrices:
            continue

        qa, pa = matrices[a]
        qn, pn = matrices[n]

        fa = features[a]
        fn = features[n]

        pair_results.append({
            "pair_id": row["pair_id"],
            "group": group,
            "delta_potassium": (
                float(row["abnormal_potassium"])
                - float(row["normal_potassium"])
            ),
            "delta_qrs_rms": (
                fa["qrs_rms"] - fn["qrs_rms"]
            ),
            "delta_post_rms": (
                fa["post_rms"] - fn["post_rms"]
            ),
            "delta_qrs_tail3": (
                fa["qrs_tail3"] - fn["qrs_tail3"]
            ),
            "delta_post_tail3": (
                fa["post_tail3"] - fn["post_tail3"]
            ),
            "delta_qrs_post_angle": (
                fa["qrs_post_angle"]
                - fn["qrs_post_angle"]
            ),
            "delta_post_window_ms": (
                fa["post_window_ms"]
                - fn["post_window_ms"]
            ),
            "qrs_matrix_rmse": float(
                np.sqrt(np.mean(
                    (qa - qn) ** 2
                ))
            ),
            "post_matrix_rmse": float(
                np.sqrt(np.mean(
                    (pa - pn) ** 2
                ))
            )
        })

pair_columns = [
    "pair_id",
    "group",
    "delta_potassium",
    "delta_qrs_rms",
    "delta_post_rms",
    "delta_qrs_tail3",
    "delta_post_tail3",
    "delta_qrs_post_angle",
    "delta_post_window_ms",
    "qrs_matrix_rmse",
    "post_matrix_rmse"
]

with RESULTS.open(
    "x", encoding="utf-8-sig", newline=""
) as f:
    writer = csv.DictWriter(
        f, fieldnames=pair_columns
    )
    writer.writeheader()
    writer.writerows(pair_results)

print("\n===== ADAPTIVE GEOMETRY V2 RESULTS =====")
print("ECGs successfully analyzed:", len(features))
print("ECGs rejected:", len(paths) - len(features))

print("\nFAILURE REASONS")
for reason, count in failures.most_common(8):
    print(reason, ":", count)

for group in (
    "HYPERKALEMIA",
    "SEVERE_HYPERKALEMIA",
    "HYPOKALEMIA"
):
    subset = [
        r for r in pair_results
        if r["group"] == group
    ]

    print("\n" + group)
    print("Original pairs:", counts[group])
    print("Pairs with valid matrices:", len(subset))

    if subset:
        print(
            "Median post-QRS RMS change:",
            round(median(
                r["delta_post_rms"] for r in subset
            ), 5)
        )
        print(
            "Median change in post-QRS window duration:",
            round(median(
                r["delta_post_window_ms"]
                for r in subset
            ), 1),
            "ms"
        )

print("\nSaved separate V2 analysis files.")
print("Original analysis files unchanged.")
print("======================================")
