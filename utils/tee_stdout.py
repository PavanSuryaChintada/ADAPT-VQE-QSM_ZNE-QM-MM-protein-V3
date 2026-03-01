"""
utils/tee_stdout.py — Tee stdout/stderr to log file
====================================================
Ensures ALL terminal output (print, logger, etc.) is saved to outputs/logs/.
"""

import sys
from datetime import datetime
from pathlib import Path


class Tee:
    """Write to both stream and a shared log file."""

    def __init__(self, stream, log_file):
        self.stream = stream
        self.log_file = log_file

    def write(self, data):
        self.stream.write(data)
        try:
            self.log_file.write(data)
            self.log_file.flush()
        except (OSError, UnicodeEncodeError):
            pass

    def flush(self):
        self.stream.flush()
        try:
            self.log_file.flush()
        except OSError:
            pass

    def close(self):
        try:
            self.log_file.close()
        except OSError:
            pass


_tee_active = None
_tee_file = None
_orig_stdout = None
_orig_stderr = None


def start_tee_session(name: str = "session") -> Path:
    """
    Start capturing all stdout/stderr to outputs/logs/{name}_{timestamp}.log.
    Call at the start of main.py, benchmark.py, etc.
    Returns the log file path.
    """
    global _tee_active, _tee_file, _orig_stdout, _orig_stderr
    if _tee_active:
        return _tee_active

    log_dir = Path("outputs/logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = log_dir / f"{name}_{ts}.log"
    _tee_file = open(log_path, "w", encoding="utf-8")

    _orig_stdout = sys.stdout
    _orig_stderr = sys.stderr
    sys.stdout = Tee(sys.stdout, _tee_file)
    sys.stderr = Tee(sys.stderr, _tee_file)
    _tee_active = log_path
    return log_path


def stop_tee_session():
    """Restore original stdout/stderr and close log file."""
    global _tee_active, _tee_file, _orig_stdout, _orig_stderr
    if _tee_active:
        try:
            if _tee_file:
                _tee_file.close()
        except OSError:
            pass
        sys.stdout = _orig_stdout or sys.__stdout__
        sys.stderr = _orig_stderr or sys.__stderr__
        _tee_active = None
        _tee_file = None
