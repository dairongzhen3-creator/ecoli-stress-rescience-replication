#!/bin/sh
set -eu
export PYTHONDONTWRITEBYTECODE=1
if [ "$#" -ne 2 ]; then echo 'Usage: sh run_source_and_figures.sh /path/to/Gervais2026_data.zip /new/output' >&2; exit 2; fi
if [ -e "$2" ]; then echo "Refusing to overwrite existing output: $2" >&2; exit 3; fi
if [ -e "$2-figures" ]; then echo "Refusing to overwrite existing figures: $2-figures" >&2; exit 3; fi
python -m unittest discover -s tests -v
python src/reviewer_run.py --author-zip "$1" --out "$2"
python src/tipping_point.py --partitions "$2/PARTITION_WEIGHT_BOUNDS/source_id_partition_counts.csv" --out "$2-tipping"
python src/replot_all_figures.py --computed-run "$2" --out "$2-figures"
