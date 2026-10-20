"""List every input of a car file with its value and source, and count what is still assumed.

    python -m cvtsim.inputs                 # the DR26 file
    python -m cvtsim.inputs my_car.yaml     # your car
    python -m cvtsim.inputs my_car.yaml --open   # only inputs still PLACEHOLDER, ASSUMED or UNVERIFIED
"""
from __future__ import annotations

import argparse
from collections import Counter

from cvtsim.io.params import DEFAULT, load_raw

TAGS = ("measured", "BRIEF", "DERIVED", "UNVERIFIED", "ASSUMED", "PLACEHOLDER")
OPEN = ("UNVERIFIED", "ASSUMED", "PLACEHOLDER")


def tag_of(source):
    """Evidence class of a source note: the strongest weakness it declares, else 'cited'."""
    s = str(source)
    for t in ("PLACEHOLDER", "ASSUMED", "UNVERIFIED"):
        if t in s:
            return t
    for t in ("DERIVED", "measured", "BRIEF"):
        if s.startswith(t):
            return t
    return "cited"


def entries(raw, prefix=""):
    """(path, value, unit, source) for every {value: ..., source: ...} entry."""
    for k, v in raw.items():
        if not isinstance(v, dict):
            continue
        path = f"{prefix}{k}"
        if "value" in v:
            yield path, v["value"], v.get("unit", ""), v.get("source", "")
        else:
            yield from entries(v, path + ".")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m cvtsim.inputs")
    ap.add_argument("car", nargs="?", default=str(DEFAULT))
    ap.add_argument("--open", action="store_true", help="show only inputs that still need data")
    a = ap.parse_args(argv)
    rows = [(p, v, u, s, tag_of(s)) for p, v, u, s in entries(load_raw(a.car))]
    count = Counter(r[4] for r in rows)
    for p, v, u, s, t in rows:
        if a.open and t not in OPEN:
            continue
        print(f"{t:<12} {p:<38} {str(v):<14} {u:<10} {s}")
    print()
    print(f"{len(rows)} inputs: " + ", ".join(f"{count[t]} {t}" for t in ("cited", *TAGS) if count[t]))
    print(f"{sum(count[t] for t in OPEN)} still need data. Run with --open to list them.")


if __name__ == "__main__":
    main()
