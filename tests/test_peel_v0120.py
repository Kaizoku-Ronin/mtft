import mpmath as mp
from mpmath import mpf, log, exp, pi
from mtft import peel, lchannels, marked_gap, coset_reps, gl2_peel, ledger_peel
from mtft import hodge_polarization as hp

def test_minf_identity():
    # K3 audit fix (BU): pin precision locally — mtft import sets dps=50,
    # but test_critical_ensemble leaves ambient dps at 20-30; these 1e-30
    # assertions must not depend on ambient precision.
    with mp.workdps(50):
        a = marked_gap.M_INF
        b = 3*log(mpf(3)/2) - log(log(mpf(3))/log(mpf(2)))
        assert abs(a-b) < mpf('1e-30')
        lam2, lam3 = log(mpf(2))/8, log(mpf(3))/27
        assert abs(log(lam2/lam3) - a) < mpf('1e-30')

def test_delta_linear_and_Rinf():
    with mp.workdps(50):  # K3 audit fix (BU): see test_minf_identity
        assert abs(marked_gap.delta_gap(3) - marked_gap.delta_gap(2)
                   - log(mpf(3)/2)) < mpf('1e-30')
        assert abs(marked_gap.R_INF - mpf('1.42507723746')) < mpf('2e-11')

def test_spectrum_ordering():
    order = [m for _, m in marked_gap.predicted_spectrum()]
    assert order == [2, 3, 5, 4, 7, 11, 9, 8, 13]

def test_bulk_expansion():
    y = mpf('0.03')
    lw = peel.w_sieve(500)
    d = peel.mu_bulk_direct(y, 500, lw)
    e = peel.mu_bulk_expansion(y, odd_terms=2)
    assert abs(d-e)/d < mpf('1e-4')

def test_skeleton_constant_is_glaisher():
    C0 = peel.skeleton_constants()['C0']
    assert abs(C0 - (1 - 12*log(mp.glaisher))) < mpf('1e-13')

def test_su5_split_identity():
    y = mpf('0.01'); lam = peel.lambda_sieve(1300)
    S, Sx, Sf = lchannels.channel_sums(5, y, lam, js=[2])
    pred = lchannels.split_formula(5, 1, S, Sx, lchannels.prime_tower(5, y))
    assert abs(Sf[1]-pred)/abs(Sf[1]) < mpf('1e-10')

def test_su13_split_identity():
    y = mpf('0.012'); lam = peel.lambda_sieve(1100)
    S, Sx, Sf = lchannels.channel_sums(13, y, lam, js=[2,4,6])
    for m in (1, 5):
        pred = lchannels.split_formula(13, m, S, Sx, lchannels.prime_tower(13, y))
        assert abs(Sf[m]-pred)/abs(Sf[m]) < mpf('1e-8')

def test_coset_accounting():
    assert coset_reps.sum_squares_check(13) and coset_reps.sum_squares_check(11)
    assert coset_reps.torus_char_count(13) == 6
    assert coset_reps.torus_char_count(11) == 5
    d = coset_reps.stage_decomposition()
    assert d['St11xSt13'] == 143 and sum(v for k,v in d.items() if k!='total') == 168

def test_conductor_certificate():
    c = gl2_peel.conductor_certificate()
    assert c['Delta'] == -1859 and c['factored'] and c['mult_11'] and c['mult_13']

def test_ap_spot():
    ap, _ = gl2_peel.ap_point_count(30)
    assert {p: ap[p] for p in (2,3,5,7,11,13,17,19,23,29)} == \
           {2:0, 3:-1, 5:-1, 7:-2, 11:-1, 13:-1, 17:-4, 19:2, 23:7, 29:-2}

def test_gl2_rank_read():
    r = gl2_peel.rank_read(mpf('8e-3'), pmax=1500)
    assert abs(r - 1) < mpf('6e-4')

def test_hodge_data():
    assert len(hp.LAMBDA_TABLE) == 11
    assert 'not-certified' in hp.LAMBDA_STATUS
    assert hp.ETA_CERT['split'] == (13, 13)

def test_ledger_addendum_verify():
    assert all(ok for _, ok in ledger_peel.verify())
