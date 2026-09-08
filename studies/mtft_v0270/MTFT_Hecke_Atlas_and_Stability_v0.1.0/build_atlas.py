"""Exact finite-field splitting atlas and independent Frobenius-degree checks."""
from collections import Counter
from fractions import Fraction
from itertools import permutations
from math import factorial, gcd, isqrt
import warnings
import sympy as sp
from sympy.polys.numberfields import galois_group
from sympy.utilities.exceptions import SymPyDeprecationWarning
from common import ROOT, X, POLYS, dump, read

warnings.filterwarnings('ignore', category=SymPyDeprecationWarning)

# Independent polynomial arithmetic: ascending coefficients, reduced modulo p.
def trim(a):
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a

def rem(a, b, p):
    a = trim([int(c) % p for c in a])
    b = trim([int(c) % p for c in b])
    inv = pow(b[-1], -1, p)
    while a != [0] and len(a) >= len(b):
        k, c = len(a)-len(b), a[-1]*inv % p
        for j, v in enumerate(b):
            a[k+j] = (a[k+j]-c*v) % p
        trim(a)
    return a

def mul(a, b, p):
    out = [0]*(len(a)+len(b)-1)
    for i, v in enumerate(a):
        for j, w in enumerate(b):
            out[i+j] = (out[i+j]+v*w) % p
    return trim(out)

def power(a, n, f, p):
    out = [1]
    while n:
        if n & 1:
            out = rem(mul(out, a, p), f, p)
        a = rem(mul(a, a, p), f, p)
        n >>= 1
    return out

def pgcd(a, b, p):
    while b != [0]:
        a, b = b, rem(a, b, p)
    return [v*pow(a[-1], -1, p) % p for v in a]

def subtract_x(a, p):
    a = a[:] + [0]*max(0, 2-len(a))
    a[1] = (a[1]-1) % p
    return trim(a)

def independent_degrees(f, p):
    n, xp, counts = len(f)-1, [0, 1], {}
    for k in range(1, n+1):
        xp = power(xp, p, f, p)
        degree = len(pgcd(f, subtract_x(xp, p), p))-1
        remaining = degree - sum(d*c for d, c in counts.items() if k % d == 0)
        assert remaining >= 0 and remaining % k == 0
        counts[k] = remaining//k
    return tuple(d for d, count in counts.items() for _ in range(count))

def irreducible_check(f, p):
    n, xp = len(f)-1, [0, 1]
    for k in range(1, n+1):
        xp = power(xp, p, f, p)
        if k < n and n % k == 0:
            assert len(pgcd(f, subtract_x(xp, p), p)) == 1
    assert rem(subtract_x(xp, p), f, p) == [0]

def partitions(n, lower=1):
    if n == 0:
        yield ()
    for first in range(lower, n+1):
        for rest in partitions(n-first, first):
            yield (first,) + rest

def density(partition):
    denom = 1
    for degree, count in Counter(partition).items():
        denom *= degree**count * factorial(count)
    return Fraction(1, denom)

def cycle_type(perm):
    seen, lengths = set(), []
    for start in range(len(perm)):
        if start in seen:
            continue
        at, length = start, 0
        while at not in seen:
            seen.add(at); length += 1; at = perm[at]
        lengths.append(length)
    return tuple(sorted(lengths))

def main():
    labels = read('inputs/prime_labels.json')
    expected_primes = [p for p in range(2, 10201) if all(p % d for d in range(2, isqrt(p)+1))]
    assert [r['p'] for r in labels] == expected_primes
    fields, cats = {}, {}
    for key, f in POLYS.items():
        group, alt = galois_group(f)
        disc = int(sp.discriminant(f))
        n = f.degree()
        assert group.order() == factorial(n) and not alt and f.is_irreducible
        # Squarefree discriminant makes Z[alpha] the full ring of integers.
        assert all(v == 1 for v in sp.factorint(disc).values())
        cats[key] = list(partitions(n))
        permutation_counts = Counter(cycle_type(p) for p in permutations(range(n)))
        assert all(Fraction(permutation_counts[c], factorial(n)) == density(c) for c in cats[key])
        fields[key] = {'polynomial': str(f.as_expr()), 'coefficients_descending': [int(v) for v in f.all_coeffs()],
                       'degree': n, 'discriminant': disc,
                       'discriminant_factorization': {str(p): int(e) for p,e in sp.factorint(disc).items()},
                       'field_discriminant_equals_polynomial_discriminant': True,
                       'galois_group': f'S{n}', 'galois_order': int(group.order()),
                       'permutation_generators': [list(g.array_form) for g in group.generators],
                       'cycle_witnesses': {}}
    assert gcd(fields['q4']['discriminant'], fields['q6']['discriminant']) == 1
    rows, checks = [], Counter()
    for label in labels:
        p = label['p']
        row = dict(label, residue_mod30=p%30, level_prime=(143 % p == 0))
        row['square_position_numerator'] = p-label['square_n']**2
        row['square_position_denominator'] = 2*label['square_n']+1
        for key, f in POLYS.items():
            _, factor_list = sp.factor_list(f, modulus=p)
            factors = [{'coefficients_descending': [int(c)%p for c in part.all_coeffs()], 'multiplicity': int(e)} for part,e in factor_list]
            factors.sort(key=lambda v:(len(v['coefficients_descending']), v['coefficients_descending']))
            product = [1]
            for part in factors:
                ff = part['coefficients_descending'][::-1]
                for _ in range(part['multiplicity']):
                    product = mul(product, ff, p)
            base = [int(c)%p for c in reversed(f.all_coeffs())]
            assert product == base
            checks['factor_products'] += 1
            pattern = tuple(sorted(len(a['coefficients_descending'])-1 for a in factors for _ in range(a['multiplicity'])))
            ramified = fields[key]['discriminant'] % p == 0
            row[key] = {'ramified':ramified, 'pattern':list(pattern), 'factors':factors}
            if not ramified:
                assert independent_degrees(base, p) == pattern
                checks['independent_degree_censuses'] += 1
                if p > 2:
                    parity = (-1)**(f.degree()-len(pattern))
                    legendre = pow(fields[key]['discriminant'] % p, (p-1)//2, p)
                    assert legendre == parity % p
                    checks['discriminant_sign_checks'] += 1
                wanted = {'irreducible':(f.degree(),), 'long_cycle':(1,f.degree()-1),
                          'transposition':(1,)*(f.degree()-2)+(2,)}
                for kind, target in wanted.items():
                    if pattern == target and kind not in fields[key]['cycle_witnesses']:
                        for part in factors:
                            irreducible_check(part['coefficients_descending'][::-1], p)
                            checks['witness_factor_irreducibility_checks'] += 1
                        fields[key]['cycle_witnesses'][kind] = {'p':p, 'pattern':list(pattern), 'factors':factors}
        row['joint_unramified'] = not(row['q4']['ramified'] or row['q6']['ramified'])
        rows.append(row)
    assert all(len(f['cycle_witnesses']) == 3 for f in fields.values())
    joint_rows = [r for r in rows if r['joint_unramified']]
    joint = []
    for a in cats['q4']:
        for b in cats['q6']:
            selected = [r for r in joint_rows if tuple(r['q4']['pattern']) == a and tuple(r['q6']['pattern']) == b]
            prob = density(a)*density(b)
            joint.append({'q4':list(a), 'q6':list(b), 'density':str(prob),
                          'limiting_percent':100*float(prob), 'count':len(selected),
                          'count_at_limiting_proportion':len(joint_rows)*float(prob),
                          'irregular_count':sum(r['irregularity_index']>0 for r in selected),
                          'regular_count':sum(r['irregularity_index']==0 for r in selected),
                          'prime_2_exception_count':sum(r['p']==2 for r in selected),
                          'first_in_sample': selected[0]['p'] if selected else None,
                          'infinitude':'proved: Chebotarev for S4 x S6',
                          'irregular_intersection_infinitude':'not established by this experiment'})
    assert len(joint) == len(cats['q4'])*len(cats['q6']) == 55
    assert sum(c['count'] for c in joint) == len(joint_rows)
    assert sum(Fraction(c['density']) for c in joint) == 1
    marginals = {}
    for key in POLYS:
        usable = [r for r in rows if not r[key]['ramified']]
        marginals[key] = [{'pattern':list(c), 'density':str(density(c)),
                           'count':sum(tuple(r[key]['pattern'])==c for r in usable)} for c in cats[key]]
    certificate = {'fields':fields, 'splitting_field_intersection':'Q', 'compositum_galois_group':'S4 x S6',
                   'compositum_degree':17280,
                   'proof_group':'Irreducible reduction implies transitivity. A (1,n-1) cycle makes a point stabilizer transitive on the other roots, hence the group is doubly transitive. Conjugating a transposition in a doubly transitive group supplies every transposition, so the group is Sn.',
                   'proof_disjointness':'The Galois closures ramify only at primes dividing their respective polynomial discriminants. These sets are disjoint. Their intersection is therefore unramified at every finite prime over Q; the Minkowski discriminant bound forces this intersection to be Q. Since both closures are Galois they are linearly disjoint, giving S4 x S6.',
                   'alternative_proof':'Any common Galois subextension gives a common quotient of S4 and S6. Their only possible nontrivial common quotient is C2. Its quadratic fields would be Q(sqrt(1957)) and Q(sqrt(194616205)), which differ. Thus the common quotient is trivial.',
                   'density_rule':'For cycle partition lambda in Sn, density = 1 / product_j(j^m_j m_j!). Joint densities multiply in the direct product. These are limiting densities, not finite-sample probabilities or short-interval coverage guarantees.',
                   'checks':dict(checks)}
    dump('field_certificate.json', certificate)
    dump('prime_atlas.json', rows)
    dump('splitting_types.json', {'marginals':marginals, 'joint':joint})
    summary = {'p_max_exclusive':10201, 'prime_count':len(rows), 'joint_unramified_count':len(joint_rows),
               'joint_excluded_primes':[r['p'] for r in rows if not r['joint_unramified']],
               'joint_types':55, 'joint_types_seen':sum(c['count']>0 for c in joint),
               'joint_types_unseen':sum(c['count']==0 for c in joint),
               'wieferich_examples':[r for r in rows if r['wieferich_base2']],
               'both_completely_split_primes':[r['p'] for r in joint_rows if r['q4']['pattern']==[1]*4 and r['q6']['pattern']==[1]*6]}
    dump('atlas_summary.json', summary)
    print({k:v for k,v in summary.items() if k!='wieferich_examples'}, flush=True)

if __name__ == '__main__':
    main()
