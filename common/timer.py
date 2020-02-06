import threading
from contextlib import ContextDecorator
import time
from typing import Any


class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""


class AggregateTime:

    def __init__(self):
        self.elapsed = 0
        self.count = 0

    def add_elapsed(self, elapsed):
        self.elapsed += elapsed
        self.count += 1

    def __str__(self):
        if self.count > 0:
            return "elapsed=%f count=%i Time per=%f" % (self.elapsed, self.count, self.elapsed / self.count)
        else:
            return ""


class Timer(ContextDecorator):
    """Time your code using a class, context manager, or decorator"""

    timers = dict()
    lock = threading.RLock()

    def __init__(self, name=None, text="Elapsed time: {:0.4f} seconds", logger=None):
        self._start_time = None
        self.name = name
        self.text = text
        self.logger = logger

        if name:
            self.lock.acquire()
            self.timers.setdefault(name, AggregateTime())
            self.lock.release()

    @classmethod
    def report_all(cls):
        cls.lock.acquire()
        all_timers = ""
        for name, agg in cls.timers.items():
            if agg.count > 0:
                all_timers += "%s: %s\n" % (name, str(agg))
        cls.lock.release()
        return all_timers

    def start(self) -> None:
        """Start a new timer"""
        if self._start_time is not None:
            raise TimerError(f"Timer is running. Use .stop() to stop it")

        self._start_time = time.perf_counter()

    def stop(self) -> float:
        """Stop the timer, and report the elapsed time"""
        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")

        # Calculate elapsed time
        elapsed_time = time.perf_counter() - self._start_time
        self._start_time = None

        # Report elapsed time
        if self.logger:
            self.logger(self.text.format(elapsed_time))
        if self.name:
            self.lock.acquire()
            self.timers[self.name].add_elapsed(elapsed_time)
            self.lock.release()

        return elapsed_time

    def __enter__(self) -> "Timer":
        """Start a new timer as a context manager"""
        self.start()
        return self

    def __exit__(self, *exc_info: Any) -> None:
        """Stop the context manager timer"""
        self.stop()
