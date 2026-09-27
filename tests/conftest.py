import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
DATA = os.path.join(ROOT, "tests", "data")
SAMPLE = os.path.join(ROOT, "Data")


def pytest_addoption(parser):
    parser.addoption("--run-slow", action="store_true", default=False, help="run slow end-to-end tests")


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: slow end-to-end test (use --run-slow)")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run-slow"):
        return
    skip = pytest.mark.skip(reason="slow; use --run-slow")
    for it in items:
        if "slow" in it.keywords:
            it.add_marker(skip)
