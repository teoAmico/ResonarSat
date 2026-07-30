#!/bin/zsh
# Create a run directory and seed its RUN.md. See runs/README.md.
#
# Exists so that recording provenance is easier than not recording it. The
# alternative -- remembering to write it afterwards -- is what produced cubes
# labelled "model unknown" and figures whose processing route had to be
# reconstructed from shell history.
#
# usage: tools/new-run.sh <scene> <suffix> "<what question this answers>"
set -e
[ $# -ge 2 ] || { echo "usage: $0 <scene> <suffix> [question]" >&2; exit 1; }
scene=$1; suffix=$2; question=${3:-"(not stated)"}
dir="runs/$scene/$(date +%Y-%m-%d)-$suffix"
mkdir -p "$dir"
{
  echo "# Run: $(date +%Y-%m-%d) $scene / $suffix"
  echo
  echo "**Question this run is meant to answer:** $question"
  echo
  echo "- git commit: \`$(git rev-parse --short HEAD 2>/dev/null || echo unknown)\`"
  echo "- started:    $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "- host:       $(uname -sm)"
  echo
  echo "## Collect"
  echo
  echo "\`\`\`"
  echo "(fill in: file, size, dwell, pulses, range bins)"
  echo "\`\`\`"
  echo
  echo "## Commands"
  echo
  echo "\`\`\`sh"
  echo "(paste each command verbatim, including every option)"
  echo "\`\`\`"
  echo
  echo "## Result"
  echo
  echo "*To be completed. A null result stays here rather than being deleted.*"
} > "$dir/RUN.md"
echo "$dir"
