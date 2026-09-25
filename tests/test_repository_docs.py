"""Navigation and generated-index contracts after moving the release history."""
import importlib.util
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_capability_index_covers_current_source():
    spec = importlib.util.spec_from_file_location(
        "build_capability_index", ROOT / "scripts/build_capability_index.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert (ROOT / "docs/CAPABILITIES.md").read_text(encoding="utf-8") == module.render()


def test_current_navigation_links_resolve():
    paths = [ROOT / "README.md", ROOT / "docs/CAPABILITIES.md",
             ROOT / "docs/SM/INTERACTION_ATLAS.md", ROOT / "docs/changelog/README.md"]
    for path in paths:
        for href in re.findall(r"\]\(([^)]+)\)", path.read_text(encoding="utf-8")):
            if href.startswith(("https://", "http://", "#")):
                continue
            assert (path.parent / href.split("#")[0]).exists(), (path, href)


def test_changelogs_are_indexed_in_one_folder():
    assert not list(ROOT.glob("CHANGELOG_*.md"))
    directory = ROOT / "docs/changelog"
    index = (directory / "README.md").read_text(encoding="utf-8")
    archived = list(directory.glob("CHANGELOG_*.md"))
    assert archived
    assert all(f"({p.name})" in index for p in archived)


def test_prebuilt_atlas_matches_source():
    from mtft.interactions import load_catalog
    from mtft.interactions.render import render_html
    assert (ROOT / "viz/sm_interaction_atlas.html").read_text(encoding="utf-8") == render_html(load_catalog())
