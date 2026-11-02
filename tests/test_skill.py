import re

import pytest

import plotnine_genomic_tracks as pgt
from plotnine_genomic_tracks._skills.install import (
    bundled_skill_dir,
    install_skill,
    main,
)

SKILL = bundled_skill_dir()


def skill_files():
    return [SKILL / "SKILL.md", *sorted((SKILL / "references").glob("*.md"))]


def test_frontmatter_and_routes():
    text = (SKILL / "SKILL.md").read_text()
    front = re.match(r"^---\n(.*?)\n---\n", text, re.S)
    assert front, "SKILL.md needs YAML frontmatter"
    fields = dict(line.split(": ", 1) for line in front.group(1).splitlines())
    assert fields["name"] == "plotnine-genomic-tracks"
    assert 100 < len(fields["description"]) <= 1024
    routed = set(re.findall(r"`(references/[\w.-]+)`", text))
    on_disk = {f"references/{p.name}" for p in (SKILL / "references").glob("*.md")}
    assert routed == on_disk, f"routed {routed} but on disk {on_disk}"


def test_mentions_only_real_api():
    import pygenomic_tracks as pygt

    pgt_names, pygt_names = set(), set()
    for path in skill_files():
        text = path.read_text()
        pgt_names |= set(re.findall(r"\bpgt\.([A-Za-z_]\w*)", text))
        pygt_names |= set(re.findall(r"\bpygt\.([A-Za-z_]\w*)", text))
        for ref in re.findall(r"`(references/[\w.-]+)`", text):
            assert (SKILL / ref).is_file(), f"{path.name} points to missing {ref}"
    assert not sorted(n for n in pgt_names if not hasattr(pgt, n))
    assert not sorted(n for n in pygt_names if not hasattr(pygt, n))


def test_documented_signatures_match():
    """Keyword arguments shown in signature blocks exist on the functions."""
    import inspect

    text = (SKILL / "references" / "tracks.md").read_text()
    for name, args in re.findall(r"pgt\.(\w+)\(([^)]*)\)\n```", text, re.S):
        params = set(inspect.signature(getattr(pgt, name)).parameters)
        shown = set(re.findall(r"(\w+)=", args))
        assert shown <= params, f"{name}: {shown - params} not parameters"


def test_install(tmp_path, capsys):
    dest = install_skill(tmp_path / "skills" / "plotnine-genomic-tracks")
    assert (dest / "SKILL.md").read_text() == (SKILL / "SKILL.md").read_text()
    assert {p.name for p in (dest / "references").glob("*.md")} == {
        p.name for p in (SKILL / "references").glob("*.md")
    }
    with pytest.raises(FileExistsError):
        install_skill(dest)
    install_skill(dest, force=True)
    assert main(["--dest", str(dest)]) == 1
    assert main(["--dest", str(dest), "--force"]) == 0
    assert main(["--print-path"]) == 0
    assert capsys.readouterr().out.strip().endswith(str(SKILL))
