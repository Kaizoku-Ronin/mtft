"""CC-34 audit: the dark-energy formula rho_Lambda / M_P^4 = delta^-6 exp(-2/alpha) (Paper 3 / Master Paper Thm 33.1), recomputed."""
import mpmath as mp
mp.mp.dps = 20
delta = mp.mpf("4.669201609102990671853"); ainv = mp.mpf("137.035999084")
pred = delta ** -6 * mp.exp(-2 * ainv)
H0 = 67.36 * 1000 / 3.0856775814913673e22; G = 6.67430e-11; c = 299792458.0; hbar = 1.054571817e-34; GeV = 1.602176634e-10
rho_L = 0.6847 * 3 * H0 ** 2 / (8 * mp.pi * G) * c ** 2 * (hbar * c) ** 3 / GeV ** 4          # GeV^4, Planck 2018
MPl = mp.sqrt(hbar * c / G) * c ** 2 / GeV; MPr = MPl / mp.sqrt(8 * mp.pi)
print("exp(-2/alpha) =", mp.nstr(mp.exp(-2 * ainv), 5), " (paper: 1.342e-119);  delta^-6 =", mp.nstr(delta ** -6, 5))
print("formula       =", mp.nstr(pred, 5), " (paper: 1.296e-123)")
print("observed      =", mp.nstr(rho_L, 4), "GeV^4 = (", mp.nstr(rho_L ** 0.25 * 1e12, 4), "meV)^4")
print("obs/M_Pl^4 (non-reduced) =", mp.nstr(rho_L / MPl ** 4, 4), " formula/obs =", mp.nstr(pred / (rho_L / MPl ** 4), 4))
print("obs/M_P^4  (reduced)     =", mp.nstr(rho_L / MPr ** 4, 4), " formula/obs =", mp.nstr(pred / (rho_L / MPr ** 4), 4))
print("sensitivity: factor 2 -> 2.01:", mp.nstr(mp.exp(-0.01 * ainv), 4), "; power 6 -> 5:", mp.nstr(delta, 4), "; H0 = 73 km/s/Mpc:", mp.nstr((73 / 67.36) ** 2, 4))
