# MDMcleaner 0.8.7 compatibility patch

This directory documents two compatibility changes applied to
MDMcleaner 0.8.7 during the thesis analysis.

The upstream MDMcleaner Python source is not redistributed here.

Instead, `prepare_mdmcleaner_patch.sh` extracts
`read_gtdb_taxonomy.py` from the local MDMcleaner 0.8.7 container and
passes it to `patch_read_gtdb_taxonomy.py`.

The transformation reproduces two changes used during the analysis:

1. SILVA database URLs were pinned to release 138.2.
2. The GTDB MD5 manifest itself was excluded from the per-file MD5
   verification loop.

For the MDMcleaner container used in the thesis:

- upstream `read_gtdb_taxonomy.py` SHA-256:
  `56fee4259214557ad6ce4003d26f3a5eb85d3374d05586655dcd78ffb3dc3696`
- patched thesis version SHA-256:
  `e9b683ce506941e89638691ac57bd66da6620d7e2a0defb8fa528df54cae69f3`

The reconstruction was validated as an exact byte-for-byte match to
the historically used patched file.

The workflow therefore preserves the analytical modification without
redistributing the complete upstream MDMcleaner source file.
