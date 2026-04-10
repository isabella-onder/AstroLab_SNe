import numpy as np
import sncosmo
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import matplotlib.pyplot as plt

SN = '2025advo'
# ── configuration ─────────────────────────────────────────────────────────────
DATA_FILE   = f"sncosmo_nonia/{SN}_input.dat"   # path to your .dat file
REDSHIFT    = 0.014946666666666667               # set your redshift here cmb
MW_EBV      =  0# 0.02078774617067834         # Milky Way E(B-V) from e.g. Schlafly+2011

# Choose one template source depending on your classification:
#   Type Ib/Ic stripped-envelope  → "nugent-sn1bc"
#   Type II-P plateau             → "nugent-sn2p"
#   Type II-L linear              → "nugent-sn2l"
#   Alternatively, SNANA templates cover IIb, Ib, Ic, II-P, II-L:
#     "snana-2004fe" (Ic), "snana-2004gq" (Ib), "snana-2004hx" (IIP), etc.
#     Full list: sncosmo.builtins.registry.retrieve_all(sncosmo.Source)
SOURCE_NAME = "nugent-sn1bc"         # change to match your type

COSMO = FlatLambdaCDM(H0=70, Om0=0.3)

# ── load data ─────────────────────────────────────────────────────────────────
data = sncosmo.read_lc(DATA_FILE)
print(f"Loaded {len(data)} photometric points across bands: {set(data['band'])}")

# ── build model with MW + host dust ───────────────────────────────────────────
dust = sncosmo.CCM89Dust()
model = sncosmo.Model(
    source=SOURCE_NAME,
    effects=[dust, dust],
    effect_names=["host", "mw"],
    effect_frames=["rest", "obs"],
)

# fix known parameters
model.set(z=REDSHIFT, mwebv=MW_EBV, mwr_v=3.1)

# ── fit ───────────────────────────────────────────────────────────────────────
# Free parameters: t0 (epoch of peak), amplitude, host dust E(B-V)
# hostr_v is typically fixed; float it only if you have dense colour coverage
result, fitted_model = sncosmo.fit_lc(
    data,
    model,
    vparam_names=["t0", "amplitude", "hostebv"],
    bounds={
        "hostebv": (0.0, 1.0),
        "hostr_v": (3.1, 3.1),   # fixed
    },
)

print("\n── fit result ───────────────────────────────────────────────────────")
print(f"  χ²  = {result.chisq:.2f}  (ndof = {result.ndof})")
print(f"  χ²ᵣ = {result.chisq / result.ndof:.2f}")
print(f"  t0  = {result.parameters[result.param_names.index('t0')]:.2f}  [MJD]")
print(f"  host E(B-V) = {result.parameters[result.param_names.index('hostebv')]:.3f}")

# ── peak apparent magnitude ───────────────────────────────────────────────────
# Evaluate in a reference band at peak — adjust 'bessellb' to suit your data
REF_BAND  = "bessellb"
t_peak    = result.parameters[result.param_names.index("t0")]
m_peak    = fitted_model.bandmag(REF_BAND, "ab", t_peak)
print(f"\n  apparent peak mag ({REF_BAND}, AB) = {m_peak:.3f}")

# ── distance modulus and absolute magnitude ───────────────────────────────────
dist_modulus = COSMO.distmod(REDSHIFT).value          # μ = 5 log10(d_L/10 pc)
M_peak       = m_peak - dist_modulus
luminosity_solar = 10 ** (-0.4 * (M_peak - 4.74))    # L_sun, using M_sun(B) ≈ 5.44
                                                       # or bolometric M_sun = 4.74

print(f"  distance modulus μ  = {dist_modulus:.3f}")
print(f"  absolute peak mag M = {M_peak:.3f}")
print(f"  luminosity          ≈ {luminosity_solar:.3e}  L_sun")

# ── plot ──────────────────────────────────────────────────────────────────────
fig = sncosmo.plot_lc(data, model=fitted_model, errors=result.errors)
plt.suptitle(
    f"{SOURCE_NAME}  |  z = {REDSHIFT}  |  M = {M_peak:.2f}  |  χ²ᵣ = {result.chisq/result.ndof:.2f}",
    fontsize=10,
)
plt.tight_layout()
plt.savefig("lightcurve_fit.pdf", dpi=150)
plt.show()
print("\nPlot saved to lightcurve_fit.pdf")