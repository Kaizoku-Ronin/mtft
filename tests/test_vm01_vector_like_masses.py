"""VM-01: the vector-like pairs are massless at tree level in our vacuum."""
from mtft.matter import vector_like_masses as V


def test_operator_selection():
    o = V.operators()
    assert not o["Sigma psi psi'"] and o["Sigma Sigma'bar psi psi'"] and o["Sigma^4"] and o["psi psi' (Dirac)"]


def test_no_tree_level_mass_and_cutoff():
    assert not V.tree_level_mass_possible()
    assert 5.4 < V.cutoff_ratio(17) < 5.8 and 4.5 < V.cutoff_ratio(19) < 4.8
