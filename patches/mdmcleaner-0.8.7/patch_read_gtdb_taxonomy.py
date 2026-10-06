#!/usr/bin/env python3
"""
Apply the compatibility modifications used with MDMcleaner 0.8.7
in the Turrialba fresh-cheese metagenomics workflow.

The upstream MDMcleaner source is not redistributed here.
Instead, this script transforms a local upstream copy of
read_gtdb_taxonomy.py.

Changes reproduced from the thesis environment:

1. Pin SILVA downloads to release 138.2.
2. Skip the GTDB MD5 manifest self-check encountered during
   database preparation.
"""

from pathlib import Path
import re
import sys


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


if len(sys.argv) != 3:
    fail(
        "usage: patch_read_gtdb_taxonomy.py "
        "<upstream_read_gtdb_taxonomy.py> "
        "<patched_output.py>"
    )

source_path = Path(sys.argv[1])
output_path = Path(sys.argv[2])

if not source_path.is_file():
    fail(f"source file not found: {source_path}")


# ------------------------------------------------------------
# Read source while retaining line endings
# ------------------------------------------------------------

text = source_path.read_text(
    encoding="utf-8"
)


# ------------------------------------------------------------
# 1. Pin SILVA to release 138.2
# ------------------------------------------------------------

replacements = (
    (
        "https://www.arb-silva.de/fileadmin/"
        "silva_databases/current",
        "https://www.arb-silva.de/archive/"
        "release_138.2",
    ),
    (
        "ftp://arb-silva.de/current",
        "https://www.arb-silva.de/archive/"
        "release_138_2",
    ),
)

for old, new in replacements:

    count = text.count(old)

    if count != 1:
        fail(
            f"expected exactly one occurrence of "
            f"{old!r}; found {count}"
        )

    text = text.replace(
        old,
        new,
        1,
    )


# ------------------------------------------------------------
# 2. Insert GTDB MD5-manifest compatibility patch
# ------------------------------------------------------------

lines = text.splitlines(
    keepends=True
)

target_code = (
    'filename = re.sub(subpattern, "", filename)'
)

target_indices = [
    i
    for i, line in enumerate(lines)
    if target_code in line
]

if len(target_indices) != 1:
    fail(
        "expected exactly one GTDB filename-normalization "
        f"line; found {len(target_indices)}"
    )

idx = target_indices[0]

target_line = lines[idx]

indent_match = re.match(
    r"^([ \t]*)",
    target_line,
)

if indent_match is None:
    fail(
        "could not determine indentation of GTDB "
        "normalization line"
    )

base_indent = indent_match.group(1)


# Determine the indentation used one level below this block
# from the nearby existing:
#
#     if fnmatch.fnmatch(filename, pattern):
#
child_indent = None

for line in lines[idx + 1:idx + 12]:

    if "if fnmatch.fnmatch(filename, pattern):" in line:

        m = re.match(
            r"^([ \t]*)",
            line,
        )

        if m:
            child_indent = m.group(1)

        break

if child_indent is None:
    fail(
        "could not infer one-level-deeper indentation "
        "from nearby code"
    )

if len(child_indent) <= len(base_indent):
    fail(
        "inferred child indentation is not deeper "
        "than the target block"
    )


# Preserve the source newline convention.
if target_line.endswith("\r\n"):
    nl = "\r\n"
else:
    nl = "\n"


patch_lines = [
    (
        base_indent
        + "# PATCH UCR/Turrialba: avoid GTDB MD5 manifest "
          "self-check"
        + nl
    ),
    (
        base_indent
        + "if fnmatch.fnmatch(filename, MD5FILEPATTERN_GTDB):"
        + nl
    ),
    (
        child_indent
        + "continue"
        + nl
    ),
]

lines[idx + 1:idx + 1] = patch_lines

patched = "".join(lines)


# ------------------------------------------------------------
# Safety checks
# ------------------------------------------------------------

checks = (
    "https://www.arb-silva.de/archive/release_138.2",
    "https://www.arb-silva.de/archive/release_138_2",
    "# PATCH UCR/Turrialba: avoid GTDB MD5 manifest self-check",
    "if fnmatch.fnmatch(filename, MD5FILEPATTERN_GTDB):",
)

for token in checks:

    if patched.count(token) != 1:
        fail(
            f"post-patch validation failed for {token!r}: "
            f"found {patched.count(token)} occurrences"
        )


output_path.parent.mkdir(
    parents=True,
    exist_ok=True,
)

output_path.write_text(
    patched,
    encoding="utf-8",
)

print(f"PATCHED_OUTPUT={output_path}")
print("SILVA_RELEASE=138.2")
print("GTDB_MD5_MANIFEST_SELF_CHECK=SKIPPED")
print("PATCH_STATUS=PASS")
