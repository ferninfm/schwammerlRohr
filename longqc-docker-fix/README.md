# LongQC Docker build fix

This is a focused fix for the `LibMambaUnsatisfiableError` in LongQC image provisioning. It moves the Python 3.9 environment into a dedicated Conda environment rather than downgrading base Python.

From the schwammerlRohr root:

```bash
cp bin/containers/dockerfiles/longqc.Dockerfile bin/containers/dockerfiles/longqc.Dockerfile.bak
cp /path/to/longqc-docker-fix/longqc.Dockerfile bin/containers/dockerfiles/longqc.Dockerfile
python3 bin/scripts/setup_containers.py --runtime docker --tools longqc
```

The build has not been run in this environment. If another dependency conflict appears, capture the full Conda solver output.

For a smaller Docker build context, add a `.dockerignore` containing `PacBio_data/`, `PacBio_output/`, `.git/`, `.snakemake/`, `*.sif`, `*.dmnd`, `__pycache__/`, and temporary backups; check whether any other Dockerfile copies additional files from the repository first.
