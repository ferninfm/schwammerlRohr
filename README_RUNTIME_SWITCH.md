# schwammerlRohr: selectable container runtime (Docker / Apptainer / Singularity)

This is a **follow-up** for the direct FASTQ-input v3 version of schwammerlRohr. It is not a standalone replacement for older versions of the repository. It modifies `bin/Snakefile`, `run_pipeline.sh`, `bin/setup.sh`, and adds `bin/scripts/docker_exec.py`, `bin/scripts/setup_containers.py`, `bin/containers/images.json`, Dockerfiles and MetaEuk/Prodigal/Tiara Apptainer definitions. Original BUSCO and mycoBinR analysis commands are retained.

## Apply it to your existing branch

From the root of your schwammerlRohr checkout, *after the external FASTQ v3 changes*:

```bash
git status --short
git apply --check /path/to/container-runtime-followup.patch
git apply /path/to/container-runtime-followup.patch
git diff --check
bash -n run_pipeline.sh bin/setup.sh
python3 -m py_compile bin/scripts/docker_exec.py bin/scripts/setup_containers.py
```

If the patch check fails, **do not overwrite your local files**. The patch is against the v3 external FASTQ code, so any different local edits may need a merge. The ZIP also contains a `changed-files/` tree for inspection, not for blindly copying over a modified checkout.

If your previous changes are staged, use `git diff --cached` for those and `git diff` for the new unstaged runtime changes. To see all changes relative to HEAD, use `git diff HEAD`.

## Use Docker (your workstation)

First, **dry-run**. This does not build images or touch your raw sequencing data:

```bash
bash run_pipeline.sh \
  --cores 16 \
  --reads '/media/fernando/data2/02_Kiruna_PacBio/Kiruna_PacBio_2025/pr_197_001.hifireads.fastq.gz' \
  --db /media/fernando/data2/nr/nr_diamond.dmnd \
  --container-runtime docker \
  --gene-prediction busco \
  --dry-run
```

For a **real run**, remove `--dry-run`. The launcher checks/builds/pulls only images required by the selected workflow. Docker images are cached by tag and reused on subsequent runs. `docker info` must work for your Linux user. You can explicitly prebuild images with `bash bin/setup.sh --runtime docker` (builds all images, including standalone predictors) or `bash bin/setup.sh --runtime docker --tools flye busco diamond report` (subset).

For standalone prediction, add `--gene-prediction standalone --metaeuk-ref /path/to/euk_ref.faa --tiara-mode run`. BUSCO and the existing report still run; the standalone DIAMOND table is an additional result, as in v3.

## Use Apptainer/Singularity (workstation or HPC)

```bash
bash run_pipeline.sh --reads /data/sample.fastq.gz \
  --db /db/nr.dmnd --container-runtime apptainer --dry-run
```

Remove `--dry-run` to let the launcher provision missing SIF files first. Use `--container-runtime singularity` if the executable is named `singularity`. Existing `.def` files remain available for the seven original tools, and new definitions are provided for MetaEuk, Prodigal and Tiara.

`--container-runtime auto` (default) prefers Apptainer, then Singularity, then Docker. `--container-runtime none` executes tools from your host PATH (you are responsible for all dependencies). `--no-container-setup` skips prebuilding images in a real run, for pre-provisioned systems.

## Image map and mounts

`bin/containers/images.json` centralizes each tool's Docker tag, SIF path, Dockerfile and Apptainer definition. Use `--container-images /path/to/custom-images.json` to use a modified copy. For published Docker images, set the `docker` reference and set `dockerfile` to `null` to pull rather than build. For Apptainer, `pull` can be a `docker://...` URI; its `apptainer` key is the destination SIF file.

In Docker mode, the repository is mounted **read/write at the same absolute path**; external FASTQ, DIAMOND and MetaEuk/Tiara reference directories are mounted **read-only**. Containers run with the current UID/GID and supplementary groups to retain output ownership. The container starts in the repository root. Paths inside tools thus remain identical to those visible from the host. Apptainer mode retains Snakemake's existing `--apptainer-args --bind` semantics.

## What is and is not validated

Local tests check the runtime CLI, read inputs with leading zeros, manifest paths, Docker command construction and mount semantics, dry-run argument generation, plus the existing v3 classification fixtures. **No real Snakemake 9.27.0 dry-run, Docker image build, Apptainer build or biological-tool execution could be performed in the authoring environment.** You must run a real dry-run and initially test image builds on your workstation; external package repositories and tool dependencies can change, and the original Singularity definition files are not all version-pinned. Dockerfiles are intended as functional equivalents of the original `.def` files, not verified bit-for-bit identical environments.

**Troubleshooting:**

- `Cannot connect to the Docker daemon`: confirm `docker info` works as your normal user.
- Docker build errors: inspect `docker build` output; try `bash bin/setup.sh --runtime docker --tools <name> --force` after resolving the issue.
- Missing SIF during dry-run: this should no longer occur solely because `report.sif` was declared as a regular Snakemake input; the `REPORT_SIF` input was removed.
- An input path not visible inside a container: check that it is declared through `--reads`, `--db`, `--metaeuk-ref`, or `--tiara-tsv`; arbitrary additional external paths require new mount configuration.
- A container executable or tool missing at runtime: validate the new Dockerfile/definition for that tool; do not assume a successful dry-run verifies container contents.

## Tests

```bash
python3 -m unittest discover -s tests -p 'test_container_runtimes.py' -v
```

The legacy workflow input interface remains unchanged: direct `--reads`, explicitly named `--sample ID=FILE`, and legacy `--samples` all work as before.
