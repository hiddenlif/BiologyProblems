"""Runs both case sets of this repository and writes their results under one output folder.

assay_math.py: receptor binding, counting beads and rare-event detection.
inference_pitfalls.py: absolute counts, technical replication and donor leakage.
"""

import json
from pathlib import Path

import assay_math
import inference_pitfalls

ROOT = Path(__file__).parent


def run(output=ROOT / "artifacts/run"):
    output = Path(output)
    return {
        "assay_math": assay_math.run(output / "assay_math"),
        "inference_pitfalls": inference_pitfalls.run(output / "inference_pitfalls"),
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
