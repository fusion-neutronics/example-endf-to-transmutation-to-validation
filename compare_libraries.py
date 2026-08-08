#!/usr/bin/env python3
"""Step 4: is one library better than another on this benchmark?

Reads two sweeps and puts them side by side, foil by foil:

    compare_libraries.py tendl-2017 tendl-2025

Each name is a folder under ``results/``, so it is whatever ``--source`` (or
``--library``) the sweep was filed under.

What "better" means here is closeness to the measurement, so the figure of merit
is ``|C/E - 1|`` and lower wins. Two things are deliberately NOT done:

* A foil missing from either sweep is excluded rather than counted as a loss.
  A library that failed to convert a foil has not thereby produced a worse
  answer for it, and comparing 73 foils against 64 would flatter whichever one
  ran to completion.
* No foil is dropped for being an outlier. A library that is wrong by 200% on
  one foil is telling you something, and hiding it behind a median would be the
  whole point of the exercise thrown away.

The measurement's own uncertainty is carried through, because most of these
foils agree within it and a difference smaller than sigma is not a result. A
foil is only called a win when the two libraries differ by more than the
measurement sigma for that foil.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import statistics

HERE = pathlib.Path(__file__).resolve().parent


def load(source: str, results_root: pathlib.Path, experiment: str) -> dict:
    """One sweep's per-foil results, keyed by foil name.

    Read from the per-foil ``fns_<case>_<experiment>.json`` files rather than
    the sweep's own ``sweep.json``, because that summary is only written when a
    sweep runs to completion. A sweep that was interrupted, or one still in
    flight, leaves the per-foil files behind and is worth comparing anyway.
    """
    root = results_root / source / "sweep"
    if not root.is_dir():
        root = results_root / source
    if not root.is_dir():
        raise SystemExit(
            f"no results for {source!r} under {results_root / source}.\n"
            f"Run: python sweep_fns.py --library {source} --source {source} ..."
        )
    rows = {}
    for path in sorted(root.glob(f"fns_*_{experiment}.json")):
        data = json.loads(path.read_text())
        case = data.get("case")
        if case is not None:
            rows[case] = data
    if not rows:
        raise SystemExit(f"no fns_*_{experiment}.json results in {root}")
    return rows


def compare(a_name, a_rows, b_name, b_rows):
    """Rows for the paired table, plus the foils only one sweep has."""
    shared = sorted(set(a_rows) & set(b_rows))
    only_a = sorted(set(a_rows) - set(b_rows))
    only_b = sorted(set(b_rows) - set(a_rows))

    rows = []
    for case in shared:
        a, b = a_rows[case], b_rows[case]
        a_err = abs(a["median_ratio"] - 1.0)
        b_err = abs(b["median_ratio"] - 1.0)
        sigma = max(a.get("median_measurement_sigma_percent", 0.0),
                    b.get("median_measurement_sigma_percent", 0.0)) / 100.0
        gap = a_err - b_err
        if abs(gap) <= sigma:
            winner = "tie"
        else:
            winner = a_name if gap < 0 else b_name
        rows.append({
            "case": case,
            "a_ratio": a["median_ratio"], "b_ratio": b["median_ratio"],
            "a_err": a_err, "b_err": b_err,
            "sigma": sigma, "winner": winner,
            })
    return rows, only_a, only_b


def table(rows, a_name, b_name):
    """The paired table, worst disagreement first."""
    width = max([len(r["case"]) for r in rows] + [4])
    head = (f"| {'foil':<{width}} | {a_name:>12} | {b_name:>12} | "
            f"{'|C/E-1| ' + a_name:>20} | {'|C/E-1| ' + b_name:>20} | closer |")
    sep = "|" + "|".join(["-" * (width + 2), "-" * 14, "-" * 14, "-" * 22, "-" * 22, "-" * 8]) + "|"
    lines = [head, sep]
    for r in sorted(rows, key=lambda r: -abs(r["a_err"] - r["b_err"])):
        closer = {"tie": "tie"}.get(r["winner"], r["winner"])
        lines.append(
            f"| {r['case']:<{width}} | {r['a_ratio']:>12.3f} | {r['b_ratio']:>12.3f} | "
            f"{r['a_err']:>20.3f} | {r['b_err']:>20.3f} | {closer:>6} |"
        )
    return "\n".join(lines)


def verdict(rows, a_name, b_name):
    """The summary that answers the question the comparison was run to ask."""
    a_wins = sum(r["winner"] == a_name for r in rows)
    b_wins = sum(r["winner"] == b_name for r in rows)
    ties = sum(r["winner"] == "tie" for r in rows)
    a_median = statistics.median(r["a_err"] for r in rows)
    b_median = statistics.median(r["b_err"] for r in rows)
    a_mean = statistics.fmean(r["a_err"] for r in rows)
    b_mean = statistics.fmean(r["b_err"] for r in rows)

    out = [
        "",
        f"{len(rows)} foils compared.",
        f"  {a_name}: median |C/E-1| {a_median:.3f}, mean {a_mean:.3f}, closer on {a_wins}",
        f"  {b_name}: median |C/E-1| {b_median:.3f}, mean {b_mean:.3f}, closer on {b_wins}",
        f"  within measurement sigma (tie): {ties}",
    ]
    better, worse = (a_name, b_name) if a_median < b_median else (b_name, a_name)
    decided = a_wins + b_wins
    if a_median == b_median or decided == 0:
        out.append("\nVerdict: nothing separates them on this benchmark.")
    elif ties >= decided:
        out.append(
            f"\nVerdict: {better} is closer on the median, but more foils are inside the "
            f"measurement uncertainty ({ties}) than are decided by it ({decided}). "
            f"That is weak evidence, not a result."
        )
    else:
        out.append(
            f"\nVerdict: {better} is closer to the measurement than {worse} on this "
            f"benchmark, on {max(a_wins, b_wins)} of {decided} decided foils."
        )
    return "\n".join(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("a", help="first results/<source> folder, e.g. tendl-2017")
    parser.add_argument("b", help="second results/<source> folder, e.g. tendl-2025")
    parser.add_argument("--experiment", default="2000exp_5min",
                        help="which campaign to compare (default: 2000exp_5min)")
    parser.add_argument("--results", type=pathlib.Path, default=HERE / "results",
                        help="where the per-source result folders live")
    parser.add_argument("--output", type=pathlib.Path, default=None,
                        help="write the table here as markdown as well as printing it")
    args = parser.parse_args()

    a_rows = load(args.a, args.results, args.experiment)
    b_rows = load(args.b, args.results, args.experiment)
    rows, only_a, only_b = compare(args.a, a_rows, args.b, b_rows)
    if not rows:
        raise SystemExit(f"no foil is in both sweeps ({args.a}: {len(a_rows)}, "
                         f"{args.b}: {len(b_rows)})")

    text = table(rows, args.a, args.b) + "\n" + verdict(rows, args.a, args.b)
    if only_a or only_b:
        text += "\n\nExcluded, present in only one sweep:"
        if only_a:
            text += f"\n  {args.a} only: {' '.join(only_a)}"
        if only_b:
            text += f"\n  {args.b} only: {' '.join(only_b)}"
    print(text)
    if args.output:
        args.output.write_text(text + "\n")
        print(f"\n-> {args.output}")


if __name__ == "__main__":
    main()
