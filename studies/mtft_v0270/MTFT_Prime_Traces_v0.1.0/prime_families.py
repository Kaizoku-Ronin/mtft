"""Prime-family definitions: theorem status is separate from finite membership.

This is a companion experiment for the MTFT workflow, not an MTFT release.
No observed count changes a family's mathematical infinitude status.
"""
from math import gcd

RESIDUES=(1,7,11,13,17,19,23,29)
MARKS=('twin_member','sophie_germain','safe','wieferich_base2')
LABELS=('Twin-prime member','Sophie Germain','Safe prime','Wieferich, base 2')
SOURCES={
 'dirichlet':{'title':'Andrew Sutherland, MIT 18.785 Lecture 17, Theorem 17.1',
   'url':'https://math.mit.edu/classes/18.785/2015fa/LectureNotes17.pdf'},
 'twin':{'title':'D. H. J. Polymath, bounded gaps retrospective',
   'url':'https://arxiv.org/abs/1409.8361'},
 'sophie_safe':{'title':'Douglas Wikstrom, On the Security of Mix-Nets and Hierarchical Group Signatures, discussion of Cunningham chains',
   'url':'https://www.csc.kth.se/utbildning/forskar/avhandlingar/doktor/2005_2006/WikstromDouglas.pdf'},
 'wieferich':{'title':'Hester Graves and Benjamin Weiss, Introduction',
   'url':'https://arxiv.org/html/2503.19144'},
 'legendre':{'title':'Marc Chamberland and Armin Straub, Weakening the Legendre Conjecture',
   'url':'https://arxiv.org/html/2602.22502v1'}
}


def catalogue():
    rows=[{'id':'all_primes','label':'All primes','infinitude_status':'PROVED_INFINITE',
           'definition':'p is a positive rational prime','justification':'Euclid theorem',
           'source_ids':['dirichlet'],'role':'total count'}]
    for residue in RESIDUES:
        assert gcd(residue,30)==1
        rows.append({'id':f'residue30_{residue}','label':f'p = {residue} mod 30',
                     'infinitude_status':'PROVED_INFINITE',
                     'definition':f'p is prime and p mod 30 equals {residue}',
                     'parameters':{'modulus':30,'residue':residue},
                     'justification':'Dirichlet theorem; gcd(residue,30)=1',
                     'source_ids':['dirichlet'],'role':'disjoint background row'})
    definitions=['p-2 or p+2 is prime; membership counts primes, not pairs',
                 '2p+1 is prime','p is odd and (p-1)/2 is prime',
                 'p is odd and 2^(p-1) = 1 mod p^2']
    for name,label,definition,source in zip(MARKS,LABELS,definitions,('twin','sophie_safe','sophie_safe','wieferich')):
        row={'id':name,'label':label,'infinitude_status':'OPEN','definition':definition,
             'role':'overlapping property','source_ids':[source]}
        if name=='wieferich_base2':row.update(parameters={'base':2},known_examples=[1093,3511])
        if name in ('sophie_germain','safe'):row['equivalent_infinitude_problem']='Sophie Germain primes <-> safe primes via p -> 2p+1'
        if name=='twin_member':row['necessary_residues_mod30']={'for_p_greater_than':7,'residues':[1,11,13,17,19,29]}
        if name=='sophie_germain':row['necessary_residues_mod30']={'for_p_greater_than':5,'residues':[11,23,29]}
        if name=='safe':row['necessary_residues_mod30']={'for_p_greater_than':11,'residues':[17,23,29]}
        rows.append(row)
    return {'schema_version':1,'status_checked_UTC':'2026-09-06','families':rows,'sources':SOURCES,
            'finite_exceptions_to_mod30_rows':[2,3,5],
            'rules':['Infinitude status belongs to the family, not an individual prime.',
                     'All family memberships in the delivered finite range are tested exactly.',
                     'Open family labels may overlap, and also belong to a proved-infinite residue row.',
                     'Infinitude of intersections and complements is not inferred.',
                     'No fit, observed member count, or absence upgrades an OPEN status.',
                     'This is a selected catalogue, not an exhaustive taxonomy of prime families.']}
