from __future__ import annotations

"""Calculate mean signal intensities for C002/C003 TIFF images and export their ratio to Excel.

This script scans all `*.tif.frames` folders in the current directory, measures the
mean pixel intensity of `*_C002T001.tif` and `*_C003T001.tif`, and writes a summary
Excel file containing the C003/C002 ratio for each folder.
"""

from pathlib import Path
import re

import numpy as np
from openpyxl import Workbook
from PIL import Image


SCRIPT_DIR = Path(__file__).resolve().parent


def natural_sort_key(text: str) -> list[object]:
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", text)]


def calculate_mean_intensity(image_path: Path) -> float:
    with Image.open(image_path) as image:
        return float(np.asarray(image, dtype=np.float64).mean())


def list_image_folders(root_dir: Path) -> list[Path]:
    image_folders = (path for path in root_dir.iterdir() if path.is_dir() and path.name.endswith(".tif.frames"))
    return sorted(image_folders, key=lambda path: natural_sort_key(path.name))


def resolve_data_root(script_dir: Path) -> Path:
    if list_image_folders(script_dir):
        return script_dir

    parent_dir = script_dir.parent
    if list_image_folders(parent_dir):
        return parent_dir

    raise FileNotFoundError("No '*.tif.frames' folders were found next to the script or in its parent directory.")


def build_summary_workbook() -> Workbook:
    workbook = Workbook()
    worksheet = workbook.active
    worksheet.title = "Summary"
    worksheet.append(
        ["folder", "c002_file", "c003_file", "c002_mean_intensity", "c003_mean_intensity", "c003_div_c002"]
    )
    worksheet.freeze_panes = "A2"
    return workbook


def append_folder_result(worksheet, folder_path: Path) -> None:
    sample_name = folder_path.name.removesuffix(".tif.frames")
    c002_path = folder_path / f"{sample_name}_C002T001.tif"
    c003_path = folder_path / f"{sample_name}_C003T001.tif"

    if not c002_path.exists() or not c003_path.exists():
        worksheet.append([sample_name, c002_path.name, c003_path.name, None, None, "missing file"])
        return

    c002_mean_intensity = calculate_mean_intensity(c002_path)
    c003_mean_intensity = calculate_mean_intensity(c003_path)
    c003_to_c002_ratio = c003_mean_intensity / c002_mean_intensity if c002_mean_intensity else None

    worksheet.append(
        [
            sample_name,
            c002_path.name,
            c003_path.name,
            c002_mean_intensity,
            c003_mean_intensity,
            c003_to_c002_ratio,
        ]
    )


def format_numeric_columns(worksheet) -> None:
    for column_name in ("D", "E", "F"):
        for cell in worksheet[column_name][1:]:
            if isinstance(cell.value, (int, float)):
                cell.number_format = "0.000000"


def main() -> None:
    data_root = resolve_data_root(SCRIPT_DIR)
    output_file = data_root / "C003_div_C002_signal_summary.xlsx"
    workbook = build_summary_workbook()
    worksheet = workbook["Summary"]

    for folder_path in list_image_folders(data_root):
        append_folder_result(worksheet, folder_path)

    format_numeric_columns(worksheet)
    workbook.save(output_file)


if __name__ == "__main__":
    main()
