"""Console + file logging setup (TDS § 33.4).

Same messages go to the screen and to results/<STATION>_<YYYYDDD>.log. WARNINGs are also collected so
that every assumption/fallback can be written to the manifest (TDS § 9.11: printed, flagged, recorded).
"""
import logging
import sys

LOGGER_NAME = "pwv_ppp"
_FMT = "%(asctime)s %(levelname)-7s %(message)s"
_DATEFMT = "%H:%M:%S"


class _WarningCollector(logging.Handler):
    """Keeps (code, message) for every WARNING so it can be recorded in the manifest."""

    def __init__(self):
        super().__init__(level=logging.WARNING)
        self.records = []

    def emit(self, record):
        msg = record.getMessage()
        code = msg.split(":", 1)[0].strip() if ":" in msg else ""
        self.records.append({"level": record.levelname, "code": code, "message": msg})


_collector = _WarningCollector()


def get():
    return logging.getLogger(LOGGER_NAME)


def setup(level="normal"):
    """Configure screen logging. level: 'quiet' | 'normal' | 'verbose'."""
    lg = get()
    lg.setLevel(logging.DEBUG)
    for h in list(lg.handlers):
        lg.removeHandler(h)
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel({"quiet": logging.WARNING, "verbose": logging.DEBUG}.get(level, logging.INFO))
    sh.setFormatter(logging.Formatter(_FMT, _DATEFMT))
    lg.addHandler(sh)
    _collector.records = []
    lg.addHandler(_collector)
    lg.propagate = False
    return lg


def add_file(path):
    """Also write every message (INFO and above, DEBUG if verbose) to a log file."""
    lg = get()
    fh = logging.FileHandler(path, mode="w", encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%Y-%m-%d %H:%M:%S"))
    lg.addHandler(fh)
    return fh


def warnings():
    """All warnings/errors emitted since setup()."""
    return list(_collector.records)


def warn(code, message):
    """Emit a coded WARNING ('CODE: message'); codes are recorded in flags/manifest."""
    get().warning("%s: %s", code, message)
