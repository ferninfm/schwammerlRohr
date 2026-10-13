#!/usr/bin/env python3
"""Provision the images defined by bin/containers/images.json.

Docker: build tagged Dockerfiles or pull a published image. Apptainer:
build local definition files or pull OCI images to .sif. Cached images are reused.
"""
import argparse
import json
import shutil
import subprocess
from pathlib import Path


def provision(repo: Path, runtime: str, manifest: dict, tools: list[str], force: bool) -> None:
    if runtime == "none":
        return
    if runtime in ("apptainer", "singularity"):
        exe = shutil.which(runtime)
        if exe is None:
            raise RuntimeError(f"'{runtime}' is not installed")
    else:
        exe = shutil.which("docker")
        if exe is None:
            raise RuntimeError("'docker' is not installed")
        subprocess.run([exe, "info", "--format", "{{.ServerVersion}}"], check=True, stdout=subprocess.DEVNULL)
    for tool in tools:
        spec = manifest[tool]
        if runtime == "docker":
            image = spec["docker"]
            if not force and subprocess.run([exe, "image", "inspect", image], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0:
                print(f"Cached Docker image: {tool} -> {image}", flush=True)
                continue
            dockerfile = spec.get("dockerfile")
            if dockerfile:
                cmd = [exe, "build", "-f", str(repo / dockerfile), "-t", image, str(repo)]
                if force:
                    cmd.insert(2, "--pull")
            else:
                cmd = [exe, "pull", image]
        else:
            dest = repo / spec["apptainer"]
            if not force and dest.is_file() and dest.stat().st_size:
                print(f"Cached SIF image: {tool} -> {dest}", flush=True)
                continue
            dest.parent.mkdir(parents=True, exist_ok=True)
            source = spec.get("pull") or (str(repo / spec["definition"]) if spec.get("definition") else None)
            if not source:
                raise ValueError(f"No Apptainer build/pull source for '{tool}'")
            cmd = [exe, "pull" if spec.get("pull") else "build", str(dest), source]
            if force and dest.exists():
                # Existing image must be replaced explicitly with --force.
                cmd.insert(2, "--force")
        print("Provisioning:", " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True, cwd=repo)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime", choices=["auto", "docker", "apptainer", "singularity", "none"], default="auto")
    parser.add_argument("--images", default="bin/containers/images.json")
    parser.add_argument("--tools", nargs="*", default=[])
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    manifest = json.loads((repo / args.images).read_text())
    runtime = args.runtime
    if runtime == "auto":
        runtime = next((tool for tool in ("apptainer", "singularity", "docker") if shutil.which(tool)), None)
        if runtime is None:
            parser.error("No supported container runtime found")
    # Preserve the original setup.sh behaviour of creating output roots.
    output_root = repo / "PacBio_output"
    for folder in ("longqc_output", "flye_output", "coverm_output", "busco_output",
                   "diamond_output", "telfinder_output", "reports", "gene_prediction"):
        (output_root / folder).mkdir(parents=True, exist_ok=True)
    tools = args.tools or list(manifest)
    bad = sorted(set(tools) - set(manifest))
    if bad:
        parser.error("Unknown tool(s): " + ", ".join(bad))
    provision(repo, runtime, manifest, tools, args.force)


if __name__ == "__main__":
    main()
