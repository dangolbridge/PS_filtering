#!/usr/bin/env python3
"""Compatibility launcher for the psfilter package CLI."""

from psfilter.cli import main


if __name__ == "__main__":
    raise SystemExit(main())
