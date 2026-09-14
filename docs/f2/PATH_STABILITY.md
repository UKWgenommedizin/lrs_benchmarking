# Path-Stability Rules

Repository reorganization must not silently break workflows or analysis scripts.

## 1. Snakemake working directory

The repository root is the canonical working directory.

```bash
cd /path/to/lrs_benchmarking
```

Workflows should build repository-local paths from the root `CWD` rather than from fragile chains such as `../../../`.

## 2. Configurable external paths

References, truth sets and other permitted external resources should accept either:

- an absolute server path; or
- a repository-relative path resolved against the project root.

Example:

```python
RAW_REFERENCE = config.get("reference", "reference/GRCh38.fa")
REFERENCE = os.path.expanduser(RAW_REFERENCE)
if not os.path.isabs(REFERENCE):
    REFERENCE = os.path.join(CWD, REFERENCE)
REFERENCE = os.path.abspath(REFERENCE)
```

Validate required paths before launching expensive jobs.

## 3. Python analysis scripts

Avoid hard-coding:

```python
Path.home() / "lrs_benchmarking"
```

because the repository may be cloned under another username, server path or container mount.

Use repository-root discovery instead:

```python
from pathlib import Path


def find_repo_root(start: Path) -> Path:
    start = start.resolve()
    for candidate in [start, *start.parents]:
        if (candidate / "CONSTITUTION.md").is_file():
            return candidate
    raise RuntimeError("Could not locate lrs_benchmarking repository root")


PROJECT_ROOT = find_repo_root(Path(__file__))
```

## 4. Before moving any existing file

Search every tracked reference first:

```bash
git grep -n "old/path/or/filename" || true
```

Then move with Git:

```bash
git mv old/path new/path
```

Update all references and re-run the search.

## 5. Compatibility entry points

For analysis scripts that are already referenced in reports or documentation, a move can leave a lightweight compatibility entry point at the previous path. This allows the repository to become organized without breaking old commands immediately.

## 6. Validation after a move

```bash
git diff --check
python3 -m py_compile path/to/script.py
```

For Snakemake workflows:

```bash
snakemake --snakefile path/to/workflow.smk --dry-run --printshellcmds
```

Do not combine large directory moves with scientific logic changes in the same commit.
