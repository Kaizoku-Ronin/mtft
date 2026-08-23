"""Gates for the RT-1B boundary arc, and smoke tests for the GP runner."""

import hashlib
import shutil
import time

import pytest

from mtft import boundary
from mtft.boundary import gates


# ---------------------------------------------------------------- provenance

def test_frozen_data_matches_provenance():
    lines = boundary.data_path("PROVENANCE.txt").read_text().splitlines()
    entries = [l.split() for l in lines if l.strip() and not l.startswith("#")]
    assert len(entries) == 5
    for sha, name in entries:
        blob = boundary.data_path(name).read_bytes()
        assert hashlib.sha1(blob).hexdigest() == sha, name


# -------------------------------------------------------------- the census

def test_boundary_dimension_formula():
    """b_N = nu_infinity + nu_2 + nu_3, recomputed from Kronecker symbols."""
    assert boundary.boundary_dim(143) == 4      # no elliptic points
    assert boundary.boundary_dim(91) == 8       # nu_3 = 4
    assert boundary.boundary_dim(65) == 8       # nu_2 = 4
    assert boundary.nu3(91) == 4 and boundary.nu2(91) == 0
    assert boundary.nu2(65) == 4 and boundary.nu3(65) == 0


def test_gate_census_consistency():
    r = gates.gate_census_consistency()
    assert r["levels"] == 22
    assert r["mismatches"] == []
    assert r["ok"]


def test_falsified_hypotheses_stay_falsified():
    """Q2 and Q4 are dead. A gate that lets them pass means the table moved."""
    r = gates.gate_hypotheses()
    assert r["computed"]["Q1"] is True
    assert r["computed"]["Q2"] is False
    assert r["computed"]["Q3"] is True
    assert r["computed"]["Q4"] is False
    assert r["ok"]


def test_q2_dies_on_the_p3_family():
    """The ten N = 3q levels give U_q = (0,0): zero total leakage."""
    p3 = [N for N, v in boundary.CENSUS.items() if v[0] == 3]
    assert len(p3) == 10
    assert all(boundary.CENSUS[N][8] == (0, 0) for N in p3)


def test_q4_dies_where_r_ell_is_two():
    big = [N for N, v in boundary.CENSUS.items() if max(v[7][1], v[8][1]) > 1]
    assert sorted(big) == [93, 111, 129, 145]


# ------------------------------------------------------- N = 143, re-derived

@pytest.mark.slow
def test_gate_h0_2k():
    r = gates.gate_h0_2k()
    assert r["dim_H0_2K"] == 36
    assert r["dim_S4"] == 40
    assert r["ok"]


@pytest.mark.slow
def test_cusp_functionals_vanish_on_h0_2k():
    r = gates.gate_cusp_functionals()
    assert r["max_value_on_H0"] == "0"
    assert r["ok"]


@pytest.mark.slow
def test_gate_leakage_exact_relations():
    """Ranks 4/2/1 with the exact vanishing support and the extra relation."""
    r = gates.gate_leakage()
    assert r["T2"]["rank"] == 4 and r["T2"]["zero_cusps"] == ()
    assert r["U11"]["rank"] == 2 and r["U11"]["zero_cusps"] == (11, 143)
    assert r["U13"]["rank"] == 1 and r["U13"]["zero_cusps"] == (13, 143)
    assert r["U13"]["extra"] == "c_1 = c_11"
    assert r["ok"]


@pytest.mark.slow
def test_vanishing_support_is_p_divisibility():
    """U_p kills exactly the cusps whose index is divisible by p."""
    r = gates.gate_leakage()
    assert all(11 % 11 == 0 for d in r["U11"]["zero_cusps"] if d % 11 == 0)
    assert set(r["U11"]["zero_cusps"]) == {d for d in (1, 11, 13, 143) if d % 11 == 0}
    assert set(r["U13"]["zero_cusps"]) == {d for d in (1, 11, 13, 143) if d % 13 == 0}


@pytest.mark.slow
def test_t2_rebuilt_from_qexpansions():
    """The shipped T_2 matrix is reproduced from the shipped S_4 basis."""
    r = gates.gate_t2_from_qexpansions()
    assert r["max_discrepancy"] == "0"
    assert r["ok"]


# ------------------------------------------------------------------- gprun

def test_find_gp_returns_path_or_none():
    gp = __import__("mtft.gprun", fromlist=["find_gp"]).find_gp()
    assert gp is None or shutil.which(gp) or gp


@pytest.mark.skipif(__import__("mtft.gprun", fromlist=["find_gp"]).find_gp() is None,
                    reason="PARI/GP not installed")
def test_job_stamps_hash_and_runs(tmp_path, monkeypatch):
    """A job must freeze its script, hash it, and stamp the output header."""
    import mtft.gprun as gprun
    monkeypatch.setattr(gprun, "RUNS", tmp_path / "runs")
    src = "print(6*7);\nquit;\n"
    job = gprun.Job(src, "smoke.gp", gprun.find_gp())
    assert job.sha256 == hashlib.sha256(src.encode()).hexdigest()
    job.start()
    for _ in range(200):
        if not job.running and job.finished:
            break
        time.sleep(0.1)
    out = job.out.read_text()
    assert "42" in out
    assert job.sha256 in out            # hash in the header AND the footer
    assert job.script.read_text() == src
    assert job.returncode == 0
    assert (job.dir / "meta.json").exists()


@pytest.mark.skipif(__import__("mtft.gprun", fromlist=["find_gp"]).find_gp() is None,
                    reason="PARI/GP not installed")
def test_partial_run_is_readable_and_marked(tmp_path, monkeypatch):
    """An interrupted job must leave readable output flagged as PARTIAL."""
    import mtft.gprun as gprun
    monkeypatch.setattr(gprun, "RUNS", tmp_path / "runs")
    job = gprun.Job("print(\"started\");\nfor(i=1,10^11, 1+1);\nquit;\n",
                    "slow.gp", gprun.find_gp())
    job.start()
    for _ in range(300):
        if job.proc is not None and "started" in job.out.read_text():
            break
        time.sleep(0.1)
    job.stop()
    # terminate -> kill can take ~10s against a tight GP loop; the PARTIAL
    # footer is only written once the process is actually reaped.
    for _ in range(400):
        if job.finished:
            break
        time.sleep(0.1)
    assert job.finished, "job never reaped after stop()"
    out = job.out.read_text()
    assert "started" in out              # flushed before termination
    assert "PARTIAL" in out              # and labelled as not a result


@pytest.mark.skipif(__import__("mtft.gprun", fromlist=["find_gp"]).find_gp() is None,
                    reason="PARI/GP not installed")
def test_stop_before_process_exists(tmp_path, monkeypatch):
    """Stop pressed before Popen returns must still cancel the job.

    Regression: stop() used to test `if self.proc` and silently do nothing
    when the subprocess had not been created yet, leaving an unreachable
    job running.
    """
    import mtft.gprun as gprun
    monkeypatch.setattr(gprun, "RUNS", tmp_path / "runs")
    job = gprun.Job("print(1);\nfor(i=1,10^12, 1+1);\nquit;\n",
                    "early.gp", gprun.find_gp())
    job.start()
    job.stop()
    for _ in range(400):
        if job.finished:
            break
        time.sleep(0.1)
    assert job.finished, "early stop was ignored"
    assert job.returncode != 0
    assert "PARTIAL" in job.out.read_text()


# ------------------------------------------------- packaging regression guard

def test_manifest_covers_every_studies_extension():
    """MANIFEST.in must ship every file type present under studies/.

    Third time this bit us. v0.16.0 dropped .txt and .json; v0.18.0 dropped
    every .gp script in studies/rt1b_2026aug -- including the frozen holdout
    runner -- because the recursive-include rule listed neither. A silently
    incomplete sdist is the failure mode that separates a claimed result
    from a reproducible one, so it gets a test rather than vigilance.
    """
    import pathlib
    root = pathlib.Path(__file__).resolve().parent.parent
    manifest = root / "MANIFEST.in"
    studies = root / "studies"
    if not manifest.exists() or not studies.exists():
        pytest.skip("running from an installed package, not a source tree")
    rule = [l for l in manifest.read_text().splitlines()
            if l.startswith("recursive-include studies")]
    assert rule, "no studies rule in MANIFEST.in"
    covered = {tok.lstrip("*") for tok in rule[0].split()[2:]}
    # __pycache__ is a build artifact, not source; everything else must ship
    present = {p.suffix for p in studies.rglob("*")
               if p.is_file() and p.suffix and "__pycache__" not in p.parts}
    missing = present - covered
    assert not missing, (
        f"file types under studies/ not shipped by MANIFEST.in: {sorted(missing)}")
