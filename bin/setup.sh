#!/usr/bin/env bash
# Provision (and reuse) only the requested per-rule images.
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec python3 "$ROOT/bin/scripts/setup_containers.py" "$@"
