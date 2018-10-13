#!/bin/bash
set -eu

# Stub example of how to run flux with library generation step.

FLUX="flux-simulator-1.2.1/bin/flux-simulator"
BUNDLE="flux-bundle"

# Run for the first time so we get .pro files
$FLUX -p $BUNDLE/"custom.par" -x

### TODO: Apply correction to .pro file after running Flux once.
### custom.pro modified
echo "Doing something to the .pro file"

# Resume with -l (library generation) and -s (sequencing) to get proper expression in FASTQ.
$FLUX -p $BUNDLE/"custom.par" -ls
