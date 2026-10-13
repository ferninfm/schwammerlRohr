#!/usr/bin/env python3
"""Run a Snakemake shell command in a per-rule Docker image.

The command is provided on stdin (quoted heredoc in Snakefile). Paths outside
workflow root are mounted read-only at their original absolute locations.
"""
import argparse
import os
from pathlib import Path
import subprocess
import sys


def docker_command(image: str, workdir: Path, mounts: list[str], script: str) -> list[str]:
    root = workdir.resolve()
    args = ["docker", "run", "--rm", "--init", "--user", f"{os.getuid()}:{os.getgid()}"]
    for gid in os.getgroups():
        if gid != os.getgid():
            args.extend(["--group-add", str(gid)])
    args += ["--mount", f"type=bind,src={root},dst={root}",
             "--workdir", str(root), "--env", "HOME=/tmp", "--env", "TMPDIR=/tmp"]
    for entry in mounts:
        directory = Path(entry).resolve(strict=True)
        if not directory.is_dir():
            raise ValueError(f"Docker mount is not a directory: {directory}")
        # Avoid duplicate/overlapping mounts when input is inside project root.
        if directory == root or root in directory.parents:
            continue
        args.extend(["--mount", f"type=bind,src={directory},dst={directory},readonly"])
    args += ["--entrypoint", "/bin/bash", image, "-e", "-u", "-o", "pipefail", "-c", script]
    return args


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", required=True)
    parser.add_argument("--mount", action="append", default=[])
    args = parser.parse_args()
    script = sys.stdin.read()
    if not script.strip():
        parser.error("empty job script")
    command = docker_command(args.image, Path.cwd(), args.mount, script)
    return subprocess.run(command, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
