"""Global cap on concurrent variant-call subprocesses across worker processes."""

from contextlib import nullcontext

_slot_sem = None


def init_worker(sem):
    global _slot_sem
    _slot_sem = sem


def set_local(sem):
    global _slot_sem
    _slot_sem = sem


def slot():
    return _slot_sem if _slot_sem is not None else nullcontext()
