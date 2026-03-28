"""Small reusable validators for public-facing configuration options."""

from __future__ import annotations

from collections.abc import Collection

from .exceptions import ConfigurationError


def ensure_choice(value: str, valid: Collection[str], *, name: str, lower: bool = True) -> str:
    """Normalize and validate one string option against an allowed set."""

    normalized = str(value).lower() if lower else str(value)
    if normalized not in valid:
        allowed = ", ".join(sorted(valid))
        raise ConfigurationError(f"{name} must be one of: {allowed}.")
    return normalized


def ensure_positive_int(value, *, name: str, minimum: int = 1) -> int:
    """Coerce one integer option and validate a lower bound."""

    integer = int(value)
    if integer < minimum:
        raise ConfigurationError(f"{name} must be at least {minimum}.")
    return integer


def ensure_probability(value, *, name: str) -> float:
    """Validate that one numeric option lies in the open interval (0, 1]."""

    numeric = float(value)
    if not 0 < numeric <= 1:
        raise ConfigurationError(f"{name} must be in the interval (0, 1].")
    return numeric


def ensure_fraction(value, *, name: str) -> float:
    """Validate that one numeric option lies in the open interval (0, 100]."""

    numeric = float(value)
    if not 0 < numeric <= 100:
        raise ConfigurationError(f"{name} must be in the interval (0, 100].")
    return numeric


def ensure_finite_float(value, *, name: str) -> float:
    """Validate that one option can be coerced to a finite float."""

    numeric = float(value)
    if not (numeric == numeric and numeric not in (float("inf"), float("-inf"))):
        raise ConfigurationError(f"{name} must be a finite number.")
    return numeric
