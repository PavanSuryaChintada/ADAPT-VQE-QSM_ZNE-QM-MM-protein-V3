"""
utils/logger.py — Custom logger (stdlib only, no loguru dependency)
====================================================================
Drop-in replacement so the code works even before loguru is installed.
Once loguru is installed, replace with: from loguru import logger
"""

import logging
import sys
import os
from datetime import datetime
from pathlib import Path


def _make_logger(name: str = "quantum_protein") -> logging.Logger:
    log = logging.getLogger(name)
    log.setLevel(logging.DEBUG)

    if log.handlers:
        return log

    # Console handler — colored
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    fmt = logging.Formatter(
        fmt="%(asctime)s | %(levelname)-8s | %(message)s",
        datefmt="%H:%M:%S"
    )
    ch.setFormatter(fmt)
    log.addHandler(ch)

    # File handler
    try:
        log_dir = Path("outputs/logs")
        log_dir.mkdir(parents=True, exist_ok=True)
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        fh = logging.FileHandler(log_dir / f"run_{ts}.log", encoding="utf-8")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging.Formatter(
            "%(asctime)s | %(levelname)s | %(message)s"
        ))
        log.addHandler(fh)
    except Exception:
        pass

    return log


class _Logger:
    """Thin wrapper giving loguru-style API over stdlib logging."""

    def __init__(self):
        self._log = _make_logger()

    def info(self, msg, *a, **k):    self._log.info(str(msg), *a, **k)
    def debug(self, msg, *a, **k):   self._log.debug(str(msg), *a, **k)
    def warning(self, msg, *a, **k): self._log.warning(str(msg), *a, **k)
    def error(self, msg, *a, **k):   self._log.error(str(msg), *a, **k)
    def success(self, msg, *a, **k): self._log.info("✓ " + str(msg), *a, **k)

    def remove(self): pass   # noop (loguru compat)
    def add(self, *a, **k): pass  # noop (loguru compat)


# Try loguru first, fall back to stdlib
try:
    from loguru import logger
except ImportError:
    logger = _Logger()
