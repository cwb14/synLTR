"""Skip markers for tests that run external programs."""
import shutil

import pytest

needs_lastz = pytest.mark.skipif(shutil.which("lastz") is None, reason="lastz (the window aligner) not on PATH")
