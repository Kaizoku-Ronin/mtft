"""Reference bookkeeping, tensor conventions, static import and export boundaries."""
import json
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from mtft.interactions import (
    load_catalog,
    new_ledger,
    select_vertices,
    summary,
    validate_catalog,
    validate_ledger,
)
from mtft.interactions.render import render_html, vertex_svg
from mtft.interactions.ufo import UFOError, import_ufo


@pytest.fixture
def catalog():
    return load_catalog()


def test_reference_counts_and_scope(catalog):
    report = summary(catalog)
    assert (report["vertices"], report["particles"], report["physical_vertices"]) == (153, 43, 71)
    assert report["provenance"]["revision"] == "abfd2c92873b3cfe4baa65c661b580e5cc07d45f"
    assert len(catalog["coverage_gaps"]) == 3
    assert "MadTeam" in catalog["provenance"]["license_text"]
    assert set(report["sectors"]) == {
        "qcd", "ew_fermion", "gauge_self", "higgs_yukawa", "higgs_gauge",
        "higgs_self", "ghost", "goldstone",
    }


def test_charge_conservation_and_fermion_flow(catalog):
    for vertex in catalog["vertices"]:
        assert sum(catalog["particles"][p]["charge"] for p in vertex["particles"]) == pytest.approx(0)
    quark, antiquark = catalog["particles"]["u"], catalog["particles"]["u__tilde__"]
    assert quark["color"] == -antiquark["color"] == 3
    assert quark["pdg_code"] == -antiquark["pdg_code"] == 2


def test_independent_known_sm_anchors(catalog):
    def with_pdgs(pdgs):
        return [v for v in catalog["vertices"] if sorted(
            catalog["particles"][p]["pdg_code"] for p in v["particles"]) == sorted(pdgs)]
    assert len(with_pdgs([25, 25, 25, 25])) == 1  # keep the quartic Higgs vertex
    gggg, = with_pdgs([21, 21, 21, 21])
    assert len(gggg["color"]) == 3 and len(gggg["lorentz"]) == 3
    assert {(c["color_index"], c["lorentz_index"]) for c in gggg["couplings"]} == {(0, 0), (1, 1), (2, 2)}
    eea, = with_pdgs([-11, 11, 22])
    assert catalog["lorentz"][eea["lorentz"][0]]["structure"] == "Gamma(3,2,1)"
    w, = with_pdgs([-12, 11, 24])
    assert "ProjM" in catalog["lorentz"][w["lorentz"][0]]["structure"]
    assert not with_pdgs([-2, 2, 25])  # source limitation is explicit


def test_filters_are_exact_and_no_data_is_shared(catalog):
    assert select_vertices(catalog, particle="11") == select_vertices(catalog, particle="e-")
    assert len(select_vertices(catalog, sector="qcd")) == 8
    assert len(select_vertices(catalog, physical_only=True)) == 71
    catalog["vertices"].clear()
    assert len(load_catalog()["vertices"]) == 153
    with pytest.raises(ValueError, match="unknown particle"):
        select_vertices(load_catalog(), particle="not-a-particle")
    with pytest.raises(ValueError, match="unknown sector"):
        select_vertices(load_catalog(), sector="gravity")


@pytest.mark.parametrize("mutation", ["duplicate", "particle", "lorentz", "slot", "coupling"])
def test_invalid_tensor_references_fail(catalog, mutation):
    v = catalog["vertices"][0]
    if mutation == "duplicate":
        catalog["vertices"].append(v)
    elif mutation == "particle":
        v["particles"][0] = "missing"
    elif mutation == "lorentz":
        v["lorentz"][0] = "FFV1"
    elif mutation == "slot":
        v["couplings"][0]["color_index"] = 90
    else:
        v["couplings"][0]["coupling"] = "missing"
    with pytest.raises(ValueError):
        validate_catalog(catalog)


def test_ledger_has_no_inferred_promotions(catalog):
    ledger = new_ledger(catalog)
    counts = validate_ledger(ledger, catalog)
    assert counts["unmapped"] == 153 and sum(counts.values()) == 153
    ledger["entries"][0]["status"] = "symmetry-compatible"
    with pytest.raises(ValueError, match="requires evidence"):
        validate_ledger(ledger, catalog)
    ledger["entries"][0]["evidence"] = ["my-audit.md#vertex-1"]
    assert validate_ledger(ledger, catalog)["symmetry-compatible"] == 1
    catalog["provenance"]["gauge"] = "different convention"
    with pytest.raises(ValueError, match="mismatch"):
        validate_ledger(ledger, catalog)


def test_ledger_requires_complete_coverage(catalog):
    ledger = new_ledger(catalog)
    ledger["entries"].pop()
    with pytest.raises(ValueError, match="missing"):
        validate_ledger(ledger, catalog)
    ledger["entries"].append(ledger["entries"][0])
    with pytest.raises(ValueError, match="duplicate"):
        validate_ledger(ledger, catalog)


def test_recorded_amplitude_gate(catalog):
    ledger = new_ledger(catalog)
    row = ledger["entries"][0]
    row.update(status="amplitude-tested", evidence=["test-only-example.md"])
    with pytest.raises(ValueError, match="at least one benchmark"):
        validate_ledger(ledger, catalog)
    row["benchmarks"] = [{
        "process": "test fixture", "energy_gev": 200, "perturbative_order": "tree",
        "observable": "test fixture", "unit": "arbitrary", "conventions": "fixture",
        "reference_artifact": "fixture-ref.json", "mtft_artifact": "fixture-mtft.json",
        "mtft_revision": "test revision", "reference_value": 1.0, "mtft_value": 1.001,
        "absolute_tolerance": 0.002, "relative_tolerance": 0,
    }]
    assert validate_ledger(ledger, catalog)["amplitude-tested"] == 1
    row["benchmarks"][0]["mtft_value"] = 1.1
    with pytest.raises(ValueError, match="fails"):
        validate_ledger(ledger, catalog)
    row["benchmarks"][0]["mtft_value"] = float("nan")
    with pytest.raises(ValueError, match="finite"):
        validate_ledger(ledger, catalog)


def test_every_vertex_svg_is_well_formed_and_labels_legs(catalog):
    for v in catalog["vertices"]:
        root = ET.fromstring(vertex_svg(v, catalog))
        labels = root.findall("{http://www.w3.org/2000/svg}text")
        assert len(labels) == len(v["particles"])
        assert labels[0].text.startswith("1:")


def test_html_embeds_roundtrippable_data_and_escapes_markup(catalog):
    catalog["scope_notes"].append('</script><script>alert("no")</script>')
    page = render_html(catalog)
    payload = page.split('<script id="atlas-data" type="application/json">')[1].split('</script>')[0]
    assert json.loads(payload)["catalog"] == catalog
    assert '<script>alert("no")' not in page
    assert '__ATLAS_DATA__' not in page and 'cdn.' not in page


@pytest.fixture
def ufo(tmp_path):
    # Independently authored tiny scalar model in the supported UFO syntax.
    contents = {
        "particles.py": "from object_library import Particle\nimport parameters as Param\nh = Particle(pdg_code=25,name='h',antiname='h',spin=1,color=1,mass=Param.MH,width=Param.ZERO,texname='h',antitexname='h',charge=0)\n",
        "vertices.py": "import particles as P\nimport lorentz as L\nimport couplings as C\nV_1 = Vertex(name='V_1',particles=[P.h,P.h,P.h],color=['1'],lorentz=[L.SSS1],couplings={(0,0):C.GC_1})\n",
        "lorentz.py": "SSS1 = Lorentz(name='SSS1',spins=[1,1,1],structure='1')\n",
        "couplings.py": "GC_1 = Coupling(name='GC_1',value='-complex(0,1)*lam',order={'QED':1})\n",
        "parameters.py": "MH = Parameter(name='MH',nature='external',type='real',value=125.,texname='m')\nZERO = Parameter(name='ZERO',nature='internal',type='real',value='0.0',texname='0')\n",
        "coupling_orders.py": "QED = CouplingOrder(name='QED',expansion_order=99,hierarchy=2)\n",
        "function_library.py": "complexconjugate = Function(name='complexconjugate',arguments=('z',),expression='z.conjugate()')\n",
    }
    for name, text in contents.items():
        (tmp_path / name).write_text(text)
    return tmp_path


def _import(path):
    return import_ufo(path, source_url="https://example.org/fixture", revision="a" * 40,
                      model_path="model", gauge="scalar fixture", restrictions="none")


def test_static_import_without_running_imports(ufo):
    with (ufo / "particles.py").open("a") as f:
        f.write("import this_module_does_not_exist\n")
    c = _import(ufo)
    assert c["vertices"][0]["sector"] == "higgs_self"
    assert len(c["provenance"]["source_files"]) == 7


@pytest.mark.parametrize("text", [
    "open('executed', 'w').write('bad')\n",
    "for i in range(3):\n    pass\n",
    "h = Particle(name=__import__('os').getcwd())\n",
    "h.name = 'changed'\n",
])
def test_dynamic_ufo_fails_without_execution(ufo, text):
    (ufo / "particles.py").write_text(text)
    with pytest.raises(UFOError):
        _import(ufo)
    assert not (ufo / "executed").exists()


def test_counterterm_model_is_not_silently_truncated(ufo):
    (ufo / "CT_vertices.py").write_text("")
    with pytest.raises(UFOError, match="counterterm"):
        _import(ufo)


def test_static_antiparticle_and_numeric_charge(ufo):
    with (ufo / "particles.py").open("a") as stream:
        stream.write("q = Particle(pdg_code=1,name='q',antiname='q~',spin=2,color=3,mass=Param.ZERO,width=Param.ZERO,texname='q',antitexname='q~',charge=-1/3)\nqbar = q.anti()\n")
    particles = _import(ufo)["particles"]
    assert particles["qbar"]["charge"] == pytest.approx(1/3)
    assert particles["qbar"]["color"] == -3
    assert particles["qbar"]["pdg_code"] == -1


def test_unresolved_mass_is_reported(ufo):
    (ufo / "parameters.py").write_text("ZERO = Parameter(name='ZERO',value='0.0')\n")
    with pytest.raises(ValueError, match="unresolved mass"):
        _import(ufo)


def test_cli_exports_offline_and_preserves_existing_evidence(tmp_path):
    def run(*args):
        return subprocess.run([sys.executable, "-m", "mtft.interactions", *map(str, args)],
                              capture_output=True, text=True, check=False)
    ledger = tmp_path / "mapping.json"
    assert run("init-ledger", "-o", ledger).returncode == 0
    before = ledger.read_bytes()
    assert run("init-ledger", "-o", ledger).returncode == 2
    assert ledger.read_bytes() == before
    assert run("validate-ledger", ledger).returncode == 0
    output = tmp_path / "index.html"
    assert run("render", "-o", output, "--ledger", ledger).returncode == 0
    assert "Standard Model interaction atlas" in output.read_text()
    assert run("show", "DOES_NOT_EXIST").returncode == 2


def test_research_links_exist_in_checkout():
    root = Path(__file__).resolve().parents[1]
    guide = json.loads((root / "src/mtft/interactions/_data/guide.json").read_text())
    for sector in guide["sectors"].values():
        for path in sector["paths"]:
            assert (root / path).is_file(), path
