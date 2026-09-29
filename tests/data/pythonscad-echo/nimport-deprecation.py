"""GUI-only: nimport() must emit DeprecationWarning (and still attempt download).

Headless builds omit nimport; this file is excluded from the echo suite when
HEADLESS=ON (see tests/CMakeLists.txt).

Uses an unreachable localhost URL so the call fails after the deprecation
path without depending on a writable user library directory or a live host.
"""

import warnings

from openscad import nimport

# Port 9 is discard; connection fails quickly without serving content.
URL = "http://127.0.0.1:9/nimport_deprecation_fixture.py"

with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always", DeprecationWarning)
    try:
        nimport(URL)
        reached_download_error = False
    except RuntimeError:
        reached_download_error = True

warned = any(issubclass(w.category, DeprecationWarning) for w in caught)
print(f"deprecation_warning: {warned}")
print(f"reached_download_error: {reached_download_error}")
