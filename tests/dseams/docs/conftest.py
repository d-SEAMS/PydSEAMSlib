# Copyright (c) d-SEAMS developers
"""Trajectories from the seams-core sources (SEAMS_CORE_SRC, else the subproject)."""

import os
from pathlib import Path

import pytest

ENGINE = Path(
    os.environ.get(
        "SEAMS_CORE_SRC",
        Path(__file__).resolve().parents[3] / "subprojects" / "seams-core",
    )
)


@pytest.fixture(scope="session")
def tip4p() -> Path:
    return ENGINE / "input" / "traj" / "exampleTraj.lammpstrj"
