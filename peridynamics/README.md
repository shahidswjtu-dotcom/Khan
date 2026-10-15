# Available peridynamics versions

This directory preserves the source available in the workspace on 2026-10-15.
Each source snapshot keeps its original relative paths and bytes. The snapshots
include historical experiments and incomplete notebook-cell fragments; they are
not all independent applications.

| Location | Purpose |
| --- | --- |
| `versions/workspace/bowen_from_scratch/` | Independent Python reconstruction; entry point `run.py`, engine `pdbeam/`, mesh builder `build.py`, tests, and earlier `fresh1` and `fresh1_atomic` versions |
| `versions/workspace/bowen_reproduction/` | Earlier `bowen_pd`, `thesis_pd`, and `printed_pd` engines; `run.py`, `run_thesis.py`, `run_printed.py`; notebooks, source cells, tests and baseline engine |
| `versions/archives/` | Separate snapshots of all nine available source archives, including original/corrected Fortran, instrumented variants and Stage 13–19 notebooks |
| `RECOVERY_MANIFEST.json` | Archive identities, original paths, included/excluded entries, duplicate identities and checksums for every recovered file |
| `builds/2026-10-15/` | Build evidence and distributable artifacts from this recovery |
| `tools/` | Reproducible integrity and Python build commands |

See [VERSIONS.md](VERSIONS.md) for the folder index and [BUILD.md](BUILD.md) to
check the sources and repeat the builds. The manifest's `version_index` contains
the detailed entry points and provenance for each archive snapshot.

The latest existing report refers to a frozen `pdbeam_v4` package and
`research/reporting/README.md`. Those sources were not present in the available
workspace or nine archives. The available `bowen_from_scratch` implementation is
an earlier baseline; it must not be relabelled as that missing version. Successful
builds and software tests do not establish agreement with the experimental
benchmark. The original project audit documents retain their scientific limits.

Generated full-run histories, large checkpoints, compiled caches and reference
papers are omitted from this source collection. Small regression fixtures are
retained where the available checks require them. The manifest lists exclusions;
continuation/analysis scripts that need omitted run data must be supplied those
original outputs separately. Historical documentation may refer to the original
workspace paths or to files omitted for this reason.

To preserve subsequent work, add a new version folder and commit it together with
updated provenance and build evidence. Keep the recovered snapshots unchanged.
