"""
This file attempts loading some environment configuration variables, and
provides a sane default otherwise.
"""

try:
    from config import NCPUS
except ImportError:
    NCPUS = 2

try:
    from config import RESULTS_PATH
except ImportError:
    RESULTS_PATH = "./results/"

try:
    from config import FPLLL_PATH
except ImportError:
    FPLLL_PATH = "./out/bin/"

try:
    from config import EMERGENCY_DUMP_PREFIX
except ImportError:
    EMERGENCY_DUMP_PREFIX = "emg_dump_"