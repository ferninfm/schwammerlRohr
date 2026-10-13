#!/usr/bin/env python3
"""Convert explicit FASTQ inputs or legacy samples to a Snakemake JSON config.

Uses JSON rather than Snakemake's --config samples=001, where YAML scalar
coercion silently converts the sample label 001 into integer 1.
"""
import argparse
import json
import re
from pathlib import Path

FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
SAMPLE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def filename_sample(path):
    name = Path(path).name
    for suffix in FASTQ_SUFFIXES:
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    else:
        raise ValueError(f"Not a FASTQ filename (expected .fastq[.gz] or .fq[.gz]): {path}")
    # Remove common *trailing* sequencing-method tags, without assuming a
    # fixed lab/project prefix such as nr_ or pr_197_.
    name = re.sub(r"(?i)\.(?:hifireads|ccs|subreads)$", "", name)
    label = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._-")
    if not label:
        raise ValueError(f"Cannot derive sample ID from {path}")
    return label


def validate_id(label):
    if not SAMPLE_RE.fullmatch(label) or ".." in label:
        raise ValueError(f"Unsafe sample ID {label!r}; use letters, digits, underscore, dot or dash")
    return label


def read_path(raw):
    path = Path(raw).expanduser().resolve(strict=True)
    if not path.is_file():
        raise ValueError(f"FASTQ is not a regular file: {raw}")
    if not str(path).lower().endswith(FASTQ_SUFFIXES):
        raise ValueError(f"FASTQ must have .fastq[.gz] or .fq[.gz] suffix: {raw}")
    if path.stat().st_size == 0:
        raise ValueError(f"FASTQ is empty: {raw}")
    return str(path)


def make_config(args):
    if args.reads and args.sample:
        raise ValueError("Do not combine --reads and --sample ID=FASTQ")
    if args.sample_id and not args.reads:
        raise ValueError("--sample-id must accompany --reads")
    if args.sample_id and len(args.sample_id) != len(args.reads):
        raise ValueError("Give one --sample-id for each FASTQ supplied to --reads")
    if (args.reads or args.sample) and args.legacy_samples_given:
        raise ValueError("Do not combine direct FASTQ inputs and --samples")

    mapping = {}
    if args.reads:
        for idx, raw in enumerate(args.reads):
            label = args.sample_id[idx] if args.sample_id else filename_sample(raw)
            label = validate_id(label)
            if label in mapping:
                raise ValueError(f"Duplicate sample ID {label!r}; specify unique --sample-id values")
            mapping[label] = read_path(raw)
    elif args.sample:
        for spec in args.sample:
            label, delim, raw = spec.partition("=")
            if not delim or not label or not raw:
                raise ValueError(f"Expected --sample ID=FASTQ, got {spec!r}")
            label = validate_id(label)
            if label in mapping:
                raise ValueError(f"Duplicate sample ID {label!r}")
            mapping[label] = read_path(raw)
    else:
        legacy = args.samples if args.samples is not None else "001,002,003"
        labels = [validate_id(s.strip()) for s in legacy.split(",") if s.strip()]
        if not labels or len(labels) != len(set(labels)):
            raise ValueError("Legacy --samples must contain unique, nonempty IDs")
        # Legacy input lookup stays in the Snakefile (PacBio_data/project+ID).
        return {"project_nr": args.project_nr if args.project_nr is not None else "nr_",
                "samples": labels, "reads": {}}
    return {"project_nr": args.project_nr if args.project_nr is not None else "",
            "samples": list(mapping), "reads": mapping}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True)
    p.add_argument("--project-nr", default=None)
    p.add_argument("--samples", default=None)
    p.add_argument("--legacy-samples-given", action="store_true")
    p.add_argument("--reads", action="append", default=[])
    p.add_argument("--sample-id", action="append", default=[])
    p.add_argument("--sample", action="append", default=[])
    args = p.parse_args()
    try:
        data = make_config(args)
    except (ValueError, OSError) as exc:
        p.error(str(exc))
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(data, indent=2) + "\n", encoding="utf8")
    for label in data["samples"]:
        print(f"  {data['project_nr']}{label} <- {data['reads'].get(label, 'PacBio_data (legacy)')}")


if __name__ == "__main__":
    main()
