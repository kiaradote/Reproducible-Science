"""Lets the tests import the workflow modules without installing anything."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scrna_workflow"))
sys.path.insert(0, str(ROOT / "tests"))


import pytest


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """A test that needs the public file in testdata/ is reported as skipped, with the fix, when the file is absent."""
    outcome = yield
    rep = outcome.get_result()
    if call.excinfo is not None and call.excinfo.errisinstance(FileNotFoundError) and "testdata" in str(call.excinfo.value):
        rep.outcome = "skipped"
        rep.longrepr = (str(item.fspath), (item.location[1] or 0), "Skipped: testdata file missing; run python fetch_testdata.py")
