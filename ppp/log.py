"""Console + file logging setup (TDS § 33.4).

Same messages go to the screen and to results/<STATION>_<YYYYDDD>.log. WARNINGs are also collected so
that every assumption/fallback can be written to the manifest (TDS § 9.11: printed, flagged, recorded).
"""
import logging
import sys
import time

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
_state = {"progress": True}


class _TqdmHandler(logging.StreamHandler):
    """Console handler that prints through tqdm.write so log lines do not break progress bars."""

    def emit(self, record):
        try:
            from tqdm import tqdm
            tqdm.write(self.format(record), file=self.stream)
        except ImportError:
            super().emit(record)


def get():
    return logging.getLogger(LOGGER_NAME)


def setup(level="normal"):
    """Configure screen logging. level: 'quiet' | 'normal' | 'verbose'."""
    lg = get()
    lg.setLevel(logging.DEBUG)
    for h in list(lg.handlers):
        lg.removeHandler(h)
    sh = _TqdmHandler(sys.stdout)
    _state["progress"] = level != "quiet"
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


# ------------------------------------------------------------------------------------------- progress bars
def progress(iterable=None, total=None, desc="", unit="it", **kw):
    """tqdm progress bar (screen only; disabled with --quiet or when tqdm is not installed)."""
    try:
        from tqdm import tqdm
    except ImportError:
        tqdm = None
    if tqdm is None or not _state["progress"]:
        return _NoBar(iterable)
    return tqdm(iterable, total=total, desc=f"    {desc}", unit=unit, leave=False, dynamic_ncols=True,
                file=sys.stderr, mininterval=0.3, **kw)


class _NoBar:
    def __init__(self, iterable=None):
        self.iterable = iterable

    def __iter__(self):
        return iter(self.iterable)

    def update(self, n=1):
        pass

    def set_postfix_str(self, *a, **k):
        pass

    def close(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


# ------------------------------------------------------------------------------------------ step checklist
def _tick(ok=True):
    sym = ("\u2714" if ok else "\u2718")
    try:
        sym.encode(sys.stdout.encoding or "ascii")
        return sym
    except (UnicodeEncodeError, LookupError):
        return "[x]" if ok else "[!]"


class Steps:
    """Numbered checklist printed at the start of a run; each step is ticked when it completes."""

    def __init__(self, items):
        self.items = list(items)          # [(key, title)]
        self.done = {}
        self.t0 = time.time()
        self.t_last = self.t0
        lg = get()
        lg.info("Processing steps:")
        for i, (_k, title) in enumerate(self.items, 1):
            lg.info("   [ ] %d. %s", i, title)

    def _idx(self, key):
        return next(i for i, (k, _t) in enumerate(self.items, 1) if k == key)

    def ok(self, key, detail=""):
        i = self._idx(key)
        now = time.time()
        self.done[key] = True
        get().info("%s %d/%d %s%s  (%.1f s)", _tick(True), i, len(self.items), self.items[i - 1][1],
                   f" - {detail}" if detail else "", now - self.t_last)
        self.t_last = now

    def skip(self, key, reason=""):
        i = self._idx(key)
        self.done[key] = None
        get().info("-  %d/%d %s skipped%s", i, len(self.items), self.items[i - 1][1], f" ({reason})" if reason else "")
        self.t_last = time.time()

    def fail(self, key=None, reason=""):
        key = key or next((k for k, _t in self.items if k not in self.done), None)
        if key is None:
            return
        i = self._idx(key)
        get().error("%s %d/%d %s failed%s", _tick(False), i, len(self.items), self.items[i - 1][1],
                    f": {reason}" if reason else "")
