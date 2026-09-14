#!/usr/bin/env python3
from __future__ import annotations

import argparse
import html
import re
import subprocess
from pathlib import Path, PurePosixPath

START = "<!-- AUTO_REPOSITORY_TREE_START -->"
END = "<!-- AUTO_REPOSITORY_TREE_END -->"

MODULE_READMES = {
    "alignment_analysis/README.md": "alignment_analysis",
    "assemblers/README.md": "assemblers",
    "assemblers/whole_genome_asm/README.md": "assemblers/whole_genome_asm",
    "assemblers/whole_genome_asm/assessment/README.md": "assemblers/whole_genome_asm/assessment",
    "assembly_analysis/README.md": "assembly_analysis",
    "variant_calling_analysis/README.md": "variant_calling_analysis",
}

# Direct files in a top-level directory are listed only when there are few of them.
# Large result-heavy folders stay compact in the root README.
ROOT_DIRECT_FILE_LIMIT = 8

# These module directories start collapsed because they tend to be large.
COLLAPSED_BY_DEFAULT = {
    "figures",
    "tables",
    "results",
    "archive",
    "assessment",
    "scripts",
    "containers",
    "config",
    "docs",
}

SECTION_RE = re.compile(
    r"(?ms)^## Repository (?:structure|explorer)\n\n"
    r"<!-- AUTO_REPOSITORY_TREE_START -->.*?"
    r"^<!-- AUTO_REPOSITORY_TREE_END -->\n*"
)

def repo_root() -> Path:
    return Path(
        subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"],
            text=True,
        ).strip()
    )

def tracked_files(root: Path) -> list[str]:
    out = subprocess.check_output(
        ["git", "-C", str(root), "ls-files"],
        text=True,
    )
    return sorted(p.strip() for p in out.splitlines() if p.strip())

def esc(value: str) -> str:
    return html.escape(value, quote=False)

def md_file_link(target: str, label: str | None = None) -> str:
    label = label or PurePosixPath(target).name
    return f"[`{esc(label)}`]({target})"

def immediate_children(files: list[str], prefix: str = ""):
    prefix = prefix.strip("/")
    needle = prefix + "/" if prefix else ""
    dirs = set()
    direct_files = []

    for path in files:
        if prefix and not path.startswith(needle):
            continue

        rel = path[len(needle):] if prefix else path
        if not rel:
            continue

        parts = rel.split("/")
        if len(parts) == 1:
            direct_files.append(parts[0])
        else:
            dirs.add(parts[0])

    return sorted(dirs), sorted(direct_files)

def count_under(files: list[str], prefix: str) -> int:
    needle = prefix.rstrip("/") + "/"
    return sum(1 for path in files if path.startswith(needle))

def render_root(files: list[str]) -> str:
    dirs, root_files = immediate_children(files)

    lines = [
        "## Repository explorer",
        "",
        START,
        "Generated from Git-tracked files. Expand only the area you need. "
        "On GitHub, press **`t`** to search tracked files by name.",
        "",
        "<details>",
        "<summary><b>Top-level files</b></summary>",
        "",
    ]

    for filename in root_files:
        lines.append(f"- {md_file_link(filename)}")

    lines += ["", "</details>", ""]

    for dirname in dirs:
        subdirs, direct_files = immediate_children(files, dirname)
        total = count_under(files, dirname)
        plural = "file" if total == 1 else "files"

        lines += [
            "<details>",
            f"<summary><b>{esc(dirname)}/</b> — {total} tracked {plural}</summary>",
            "",
        ]

        readme = f"{dirname}/README.md"
        if readme in files:
            lines.append(f"- Documentation: {md_file_link(readme, 'README.md')}")

        for subdir in subdirs:
            n = count_under(files, f"{dirname}/{subdir}")
            p = "file" if n == 1 else "files"
            lines.append(f"- `{esc(subdir)}/` — {n} {p}")

        visible_direct = [
            filename
            for filename in direct_files
            if f"{dirname}/{filename}" != readme
        ]

        if visible_direct:
            if len(visible_direct) <= ROOT_DIRECT_FILE_LIMIT:
                for filename in visible_direct:
                    lines.append(
                        f"- {md_file_link(f'{dirname}/{filename}')}"
                    )
            else:
                lines.append(
                    f"- `{len(visible_direct)} direct files` — "
                    "use GitHub's **`t`** file finder or the module README to locate a specific file"
                )

        lines += ["", "</details>", ""]

    lines.append(END)
    return "\n".join(lines)

def build_tree(files: list[str], module_root: str):
    prefix = module_root.rstrip("/") + "/"
    tree = {"dirs": {}, "files": []}

    for path in files:
        if not path.startswith(prefix):
            continue

        rel = path[len(prefix):]
        if not rel:
            continue

        node = tree
        parts = rel.split("/")

        for part in parts[:-1]:
            node = node["dirs"].setdefault(
                part,
                {"dirs": {}, "files": []},
            )

        node["files"].append(parts[-1])

    return tree

def recursive_count(node) -> int:
    return len(node["files"]) + sum(
        recursive_count(child)
        for child in node["dirs"].values()
    )

def render_node(node, rel_dir: str = "") -> list[str]:
    lines = []

    for filename in sorted(node["files"]):
        target = f"{rel_dir}/{filename}" if rel_dir else filename
        lines.append(f"- {md_file_link(target)}")

    for dirname in sorted(node["dirs"]):
        child = node["dirs"][dirname]
        child_rel = f"{rel_dir}/{dirname}" if rel_dir else dirname
        n = recursive_count(child)
        plural = "file" if n == 1 else "files"
        open_attr = "" if dirname in COLLAPSED_BY_DEFAULT else " open"

        lines += [
            "",
            f"<details{open_attr}>",
            f"<summary><b>{esc(dirname)}/</b> — {n} {plural}</summary>",
            "",
        ]

        lines.extend(render_node(child, child_rel) or ["_No tracked files._"])
        lines += ["", "</details>"]

    return lines

def render_module(files: list[str], module_root: str) -> str:
    tree = build_tree(files, module_root)

    lines = [
        "## Repository explorer",
        "",
        START,
        "Generated from Git-tracked files. Expand only the directory you need. "
        "On GitHub, press **`t`** for fast filename search.",
        "",
    ]

    lines.extend(render_node(tree))
    lines += ["", END]

    return "\n".join(lines)

def replace_or_insert(text: str, generated: str) -> str:
    replacement = generated.rstrip() + "\n\n"

    if SECTION_RE.search(text):
        return SECTION_RE.sub(replacement, text, count=1)

    lines = text.splitlines()
    insert_at = None

    for index, line in enumerate(lines):
        if index > 0 and line.startswith("## "):
            insert_at = index
            break

    if insert_at is None:
        return text.rstrip() + "\n\n" + generated.rstrip() + "\n"

    before = "\n".join(lines[:insert_at]).rstrip()
    after = "\n".join(lines[insert_at:]).lstrip()

    return before + "\n\n" + generated.rstrip() + "\n\n" + after + "\n"

def update_file(path: Path, generated: str, check: bool) -> bool:
    if not path.exists():
        return False

    old = path.read_text(encoding="utf-8")
    new = replace_or_insert(old, generated)

    if new == old:
        return False

    if check:
        print(f"OUTDATED: {path}")
        return True

    path.write_text(new, encoding="utf-8")
    print(f"Updated: {path}")
    return True

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Update collapsible repository navigation in README files."
    )
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    root = repo_root()
    files = tracked_files(root)
    changed = False

    changed |= update_file(
        root / "README.md",
        render_root(files),
        args.check,
    )

    for readme, module_root in MODULE_READMES.items():
        changed |= update_file(
            root / readme,
            render_module(files, module_root),
            args.check,
        )

    if args.check:
        if changed:
            print("Repository navigation is out of date.")
            return 1

        print("Repository navigation is up to date.")
        return 0

    if not changed:
        print("Repository navigation is already up to date.")

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
