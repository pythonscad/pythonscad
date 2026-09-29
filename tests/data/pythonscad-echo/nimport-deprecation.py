"""GUI-only: nimport() must emit DeprecationWarning (and still attempt download).

Headless builds omit nimport; this file is excluded from the echo suite when
HEADLESS=ON (see tests/CMakeLists.txt).

Uses a non-HTTP(S) URL so curl_download rejects it immediately (protocols are
restricted to http/https), after the deprecation path runs, without depending
on a writable user library directory or network reachability.
"""

import warnings

from openscad import nimport

# curl_download only allows http/https; ftp is rejected before any connect.
URL = "ftp://example.com/nimport_deprecation_fixture.py"

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
