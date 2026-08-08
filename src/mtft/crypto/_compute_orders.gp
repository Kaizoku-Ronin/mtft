\\ MTFT-AL-Commit v0.1: Jacobian Order Analysis for X_0(143)
\\ Compute |A_i(F_q)| via |A_f(F_q)| = N_{K_f/Q}(q+1 - a_q)

default(parisize, 4*10^9);
default(realprecision, 38);

print("══════════════════════════════════════════════════════════════");
print(" MTFT-AL-Commit v0.1 — Jacobian Order Analysis");
print(" Level N=143, three Galois orbits A_1 (dim 1), A_2 (dim 4), A_3 (dim 6)");
print("══════════════════════════════════════════════════════════════");

N = 143;
mf = mfinit([N, 2], 0);
basis = mfeigenbasis(mf);
print("Newforms: ", #basis);

orbit_dims = vector(#basis);
field_polys = vector(#basis);
{
  for(i = 1, #basis,
    a2 = mfcoef(basis[i], 2);
    if(type(a2) == "t_POLMOD",
      orbit_dims[i] = poldegree(component(a2, 1));
      field_polys[i] = component(a2, 1)
    ,
      orbit_dims[i] = 1;
      field_polys[i] = 0
    );
    print("  Orbit ", i, ": dim = ", orbit_dims[i]);
  );
}
print("  Sum: ", sum(i = 1, #basis, orbit_dims[i]));

\\ |A_f(F_q)| = N_{K/Q}(q+1 - a_q) via resultant
norm_eval(q, aq) = if(type(aq) == "t_POLMOD", abs(polresultant(component(aq,1), (q+1) - component(aq,2))), abs(q + 1 - aq));

primes_to_test = [2,3,5,7,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113];

n_primes = #primes_to_test;
orders = matrix(n_primes, 3);

print("");
print("    q  |  |A_1|  |       |A_2|       |             |A_3|             ");
print("-------+---------+--------------------+----------------------------------");
{
  for(i = 1, n_primes,
    q = primes_to_test[i];
    for(j = 1, #basis,
      orders[i, j] = norm_eval(q, mfcoef(basis[j], q));
    );
    printf("  %4d | %7s | %18s | %s\n", q, Str(orders[i,1]), Str(orders[i,2]), Str(orders[i,3]));
  );
}

\\ Sanity check vs ellap on 143a1
print("");
print("Sanity (A_1 vs 143a1):");
E = ellinit([0, -1, 1, -1, -2]);
{
  ok = 1;
  for(i = 1, n_primes,
    q = primes_to_test[i];
    nq_curve = q + 1 - ellap(E, q);
    if(nq_curve != orders[i, 1], ok = 0; print("  MISMATCH q=", q, ": curve=", nq_curve, " vs A_1=", orders[i, 1]));
  );
  if(ok, print("  All ", n_primes, " primes match."));
}

print("");
print("══════════════════════════════════════════════════════════════");
print(" Factorizations & cryptographic assessment");
print("══════════════════════════════════════════════════════════════");

\\ CSV output
write("/home/claude/work/jacobian_data.csv", "q,orbit,dim,order,largest_prime,bits,embedding_degree");
{
  for(i = 1, n_primes,
    q = primes_to_test[i];
    for(j = 1, 3,
      n = orders[i, j];
      if(n > 1,
        fct = factor(n);
        nrows = matsize(fct)[1];
        largest = fct[nrows, 1];
        bits = log(largest)/log(2);
        emb_deg = if(largest > q && gcd(q, largest) == 1, znorder(Mod(q, largest)), 0);
        write("/home/claude/work/jacobian_data.csv", q, ",", j, ",", orbit_dims[j], ",", n, ",", largest, ",", bits, ",", emb_deg);
      );
    );
  );
}

print("CSV written to /home/claude/work/jacobian_data.csv");

\\ Print key cryptographic stats
print("");
print("Top candidates (by largest prime factor in A_3):");
print("");
print("    q  | A_3 largest prime |  bits  | emb deg");
print("-------+---------------------+--------+--------");
{
  for(i = 1, n_primes,
    q = primes_to_test[i];
    n3 = orders[i, 3];
    fct = factor(n3);
    largest = fct[matsize(fct)[1], 1];
    bits = log(largest)/log(2);
    emb_deg = if(largest > q && gcd(q, largest) == 1, znorder(Mod(q, largest)), 0);
    printf("  %4d | %19s | %6.2f | %d\n", q, Str(largest), bits, emb_deg);
  );
}

print("");
print("Detailed factorizations at q = 47, 71, 97, 113:");
{
  for(idx = 1, n_primes,
    q = primes_to_test[idx];
    if(q == 47 || q == 71 || q == 97 || q == 113,
      print("");
      print("--- q = ", q, " ---");
      for(j = 1, 3,
        n = orders[idx, j];
        print("|A_", j, "(F_", q, ")| = ", n, "  factorization: ", factor(n));
      );
    );
  );
}

quit;
