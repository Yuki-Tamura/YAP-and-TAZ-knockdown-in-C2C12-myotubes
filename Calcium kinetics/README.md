# Ca2+ Imaging Analysis for C2C12 Myotubes

This repository contains scripts used for:

1. creating aligned 2x2 comparison videos from time-lapse TIFF image sequences
2. calculating Ca2+ kinetics from fluorescence time-series data stored in Excel workbooks

The analysis was performed on C2C12 myotubes under electrical pulse stimulation and was used to compare `siScr`, `siYAP`, `siTAZ`, and `siYAPsiTAZ` conditions.

## Files

- `create_aligned_calcium_response_video.py`
  - Creates an aligned 2x2 MP4 video from four TIFF frame folders.
- `calculate_calcium_kinetics_from_excel.py`
  - Extracts time-series traces from the Excel analysis workbook and calculates per-replicate Ca2+ kinetics.
- `Analysis_20250326_210022.xlsx`
  - Excel workbook containing the fluorescence time-series data used for kinetics analysis.

## Requirements

- Python 3.13.1
- `numpy` 2.2.0
- `imageio` 2.37.3
- `Pillow` 12.2.0
- `openpyxl` 3.1.5
- FFmpeg 7.1

Install Python packages with:

```bash
python3 -m pip install numpy imageio pillow openpyxl
```

## Input Data Layout

Each movie condition is expected to be stored in a folder named:

```text
<sample_name>.tif.frames
```

For example:

```text
siScr-1.tif.frames
siYAP-1.tif.frames
siTAZ-1.tif.frames
siYAPsiTAZ-1.tif.frames
```

Each folder should contain sequential TIFF frames named like:

```text
siScr-1_T001.tif
siScr-1_T002.tif
...
siScr-1_T500.tif
```

## 1. Create an Aligned Comparison Video

`create_aligned_calcium_response_video.py` creates a 2x2 composite video using the original TIFF images.

### What it does

- loads four image sequences
- detects repeated peaks from the mean green-channel signal
- applies user-specified starting frames
- creates a labeled 2x2 comparison movie

### Example

```bash
python3 create_aligned_calcium_response_video.py \
  --panels siScr-1 siYAP-1 siTAZ-1 siYAPsiTAZ-1 \
  --starts 20 19 21 19 \
  --output aligned_peak_quad_full.mp4 \
  --title "Ca2+ Responses of C2C12 Myotubes to Electrical Pulse Stimulation"
```

### Arguments

- `--panels`
  - Four sample labels.
- `--starts`
  - Four starting frame numbers, 1-based, corresponding to the samples in `--panels`.
- `--output`
  - Output MP4 file name.
- `--title`
  - Title displayed at the top of the movie.

## 2. Calculate Ca2+ Kinetics from Excel

`calculate_calcium_kinetics_from_excel.py` reads the normalized time-series traces from the Excel workbook and writes replicate-level metrics to CSV.

### What it calculates

For each replicate:

- `rise_slope_mean`
- `decay_slope_mean`
- `time_to_peak_mean_s`
- `decay_tau_mean_s`
- `n_events`

### Definitions

- `rise_slope_mean`
  - Mean slope from the preceding trough to the peak.
- `decay_slope_mean`
  - Mean slope from the peak to the following trough.
- `time_to_peak_mean_s`
  - Mean time from the preceding trough to the peak.
- `decay_tau_mean_s`
  - Mean apparent decay time constant, defined here as the time for the signal to fall to `1/e` of the peak-to-trough amplitude during the decay phase.

### Example

```bash
python3 calculate_calcium_kinetics_from_excel.py
```

This writes:

```text
ca_kinetics_per_replicate.csv
```

## Output CSV

The output CSV contains the following columns:

```text
replicate,group,n_events,rise_slope_mean,decay_slope_mean,time_to_peak_mean_s,decay_tau_mean_s
```

## Notes

- The video script uses the original TIFF frames directly to minimize image-quality loss.
- The kinetics analysis currently uses the `Sheet1` time-series data in `Analysis_20250326_210022.xlsx`.
- Peak and trough detection is based on local extrema with a minimum gap of 5 frames.
- `decay_tau_mean_s` is an apparent decay constant derived from the observed trace and is not a full exponential-fit parameter.
