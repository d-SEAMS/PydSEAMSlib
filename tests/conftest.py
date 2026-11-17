import os
from pathlib import Path

import pytest

_SEAMS_KEYS = (
    "SEAMS_FRAME",
    "SEAMS_CUTOFF",
    "SEAMS_K",
    "SEAMS_GRAPH",
    "SEAMS_CONFIG",
)


@pytest.fixture(scope="session")
def example_traj():
    return Path(__file__).resolve().parent / "data" / "exampleTraj.lammpstrj"


@pytest.fixture(autouse=True)
def _isolate_seams_config():
    from pydseams import config

    saved = {key: os.environ[key] for key in _SEAMS_KEYS if key in os.environ}
    for key in _SEAMS_KEYS:
        os.environ.pop(key, None)
    config.reset()
    yield
    config.reset()
    for key in _SEAMS_KEYS:
        os.environ.pop(key, None)
    os.environ.update(saved)
