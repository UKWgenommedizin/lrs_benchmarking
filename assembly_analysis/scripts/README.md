# Assembly Analysis Scripts

Use two categories:

```text
metrics/   # metric extraction / table aggregation
plots/     # figure generation
```

Scripts should discover the repository root robustly rather than assuming `/home/<user>/lrs_benchmarking`.

Recommended Python pattern:

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

Then build paths from `PROJECT_ROOT`.
