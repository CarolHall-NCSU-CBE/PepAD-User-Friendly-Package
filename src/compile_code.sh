#!/bin/bash
#
# Source layout:
# src/
# |-- compile_code.sh
# |-- Dockerfile
# |-- main_v1.42-2.f90
# |-- lib/
# `-- man/
#     `-- man1/
#         `-- PepAD.1
#
# Default installation layout:
# src/
# `-- PepAD/
#     |-- PepAD
#     |-- lib/
#     `-- man/
#         `-- man1/
#             `-- PepAD.1

set -euo pipefail

# Locate the directory containing this script.
SOURCE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

# Use src/PepAD by default, or accept a user-defined installation directory.
OUTPUT_DIR="${1:-$SOURCE_DIR/PepAD}"

SOURCE_FILE="$SOURCE_DIR/main_v1.42-2.f90"
LIB_DIR="$SOURCE_DIR/lib"
MAN_FILE="$SOURCE_DIR/man/man1/PepAD.1"

# Check all required source files before changing an existing installation.
if [[ ! -f "$SOURCE_FILE" ]]; then
    echo "ERROR: Fortran source file not found: $SOURCE_FILE" >&2
    exit 1
fi

if [[ ! -d "$LIB_DIR" ]]; then
    echo "ERROR: PepAD runtime library not found: $LIB_DIR" >&2
    exit 1
fi

if [[ ! -f "$MAN_FILE" ]]; then
    echo "ERROR: PepAD manual page not found: $MAN_FILE" >&2
    exit 1
fi

if [[ -z "$OUTPUT_DIR" || "$OUTPUT_DIR" == "/" || "$OUTPUT_DIR" == "$SOURCE_DIR" ]]; then
    echo "ERROR: Unsafe installation directory: $OUTPUT_DIR" >&2
    exit 1
fi

# Load the Intel compiler environment when the module command is available.
if command -v module >/dev/null 2>&1; then
    module load PrgEnv-intel/2024.1.0
fi

if ! command -v ifx >/dev/null 2>&1; then
    echo "ERROR: Intel ifx compiler is not available." >&2
    exit 1
fi

# Preserve an existing installation before creating the new one.
if [[ -e "$OUTPUT_DIR" ]]; then
    BACKUP_DIR="${OUTPUT_DIR}_backup_$(date +%Y%m%d_%H%M%S)"
    mv -- "$OUTPUT_DIR" "$BACKUP_DIR"
    echo "Previous installation moved to: $BACKUP_DIR"
fi

mkdir -p "$OUTPUT_DIR/man/man1"

# Compile the current PepAD source code.
ifx -O2 -o "$OUTPUT_DIR/PepAD" "$SOURCE_FILE"

# Keep the runtime library beside the PepAD executable.
cp -a "$LIB_DIR" "$OUTPUT_DIR/lib"

# Install the manual page inside the PepAD installation directory.
install -m 644 "$MAN_FILE" "$OUTPUT_DIR/man/man1/PepAD.1"

echo
echo "PepAD installation completed:"
echo "  Executable: $OUTPUT_DIR/PepAD"
echo "  Library:    $OUTPUT_DIR/lib"
echo "  Manual:     $OUTPUT_DIR/man/man1/PepAD.1"
echo
echo "Use PepAD in the current shell with:"
echo "  export PATH=\"$OUTPUT_DIR:\$PATH\""
echo "  export MANPATH=\"$OUTPUT_DIR/man:\${MANPATH:-}\""
echo
echo "Test the installation with:"
echo "  PepAD --help"
echo "  man PepAD"
