"""
prodfunc.py - Python implementation of the productionRate function.

Loads reaction data from YAML files (one per reaction) and replicates the
C++ productionRate(wdot, sc, T, Te, EN, enerExch) function.

Usage example
-------------
    from prodfunc import ReactionDB
    import numpy as np

    db = ReactionDB()          # loads from tools/reactions/ by default
    wdot = np.zeros(44)
    sc   = np.zeros(44)
    ener_exch = np.zeros(44)
    # ... populate sc ...
    db.production_rate(wdot, sc, T=300.0, Te=300.0, EN=0.0, ener_exch=ener_exch)
"""

import math
import os
import glob

import yaml
import numpy as np

# ---------------------------------------------------------------------------
# Species names (indices 0-43)
# ---------------------------------------------------------------------------

SPECIES_NAMES = [
    'E',      # 0  - electron
    'H',      # 1
    'Hn',     # 2  - H-
    'H2',     # 3
    'H2v',    # 4  - H2 vibrational
    'H2p',    # 5  - H2+
    'H3p',    # 6  - H3+
    'O',      # 7
    'On',     # 8  - O-
    'O2',     # 9
    'O2n',    # 10 - O2-
    'CO',     # 11
    'CO2',    # 12
    'CO2v1',  # 13
    'CO2v2',  # 14
    'CO2v3',  # 15
    'CO2v4',  # 16
    'CO2p',   # 17 - CO2+
    'C2O3p',  # 18 - C2O3+
    'C2O4p',  # 19 - C2O4+
    'CO3n',   # 20 - CO3-
    'C',      # 21
    'OH',     # 22
    'OHp',    # 23 - OH+
    'H2O',    # 24
    'H2Op',   # 25 - H2O+
    'H3Op',   # 26 - H3O+
    'HCO',    # 27
    'CHp',    # 28 - CH+
    'CH2',    # 29
    'CH2p',   # 30 - CH2+
    'CH3',    # 31
    'CH3p',   # 32 - CH3+
    'CH4',    # 33
    'CH2O',   # 34
    'CH3O',   # 35
    'CH3OH',  # 36
    'CH2OH',  # 37
    'AR',     # 38
    'ARe',    # 39 - Ar*
    'ARp',    # 40 - Ar+
    'ARHp',   # 41 - ArH+
    'AR2e',   # 42 - Ar2*
    'AR2p',   # 43 - Ar2+
]

NUM_SPECIES = 44
NA = 6.02214085774e23   # Avogadro's number


# ---------------------------------------------------------------------------
# Stub functions – replace with actual implementations when available
# ---------------------------------------------------------------------------

def gibbs(g_RT, T):
    """
    Stub: compute species Gibbs free energy g_RT[i] = G_i / (R*T).

    This must be replaced with the actual thermodynamic implementation.
    The function should populate the length-44 array g_RT in-place.
    """
    for i in range(NUM_SPECIES):
        g_RT[i] = 0.0


def comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx,
                   ener_exch, Ue, T, Te):
    """
    Stub: compute electron energy exchange contribution from a reaction.

    Arguments mirror the C++ signature::

        comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx,
                       enerExch, Ue, tc[1], Te)

    Replace with the actual implementation.
    """
    pass


# ---------------------------------------------------------------------------
# ReactionDB
# ---------------------------------------------------------------------------

class ReactionDB:
    """
    Load all reaction YAML files and provide a production_rate() method that
    replicates productionRate() from prodfunc.cpp.
    """

    def __init__(self, reactions_dir=None):
        if reactions_dir is None:
            reactions_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         'reactions')
        self.reactions = []
        for path in sorted(glob.glob(os.path.join(reactions_dir, 'reaction_*.yaml'))):
            with open(path) as fh:
                self.reactions.append(yaml.safe_load(fh))

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def production_rate(self, wdot, sc, T, Te, EN, ener_exch):
        """
        Compute chemical production rates (mol / m³ / s) for all 44 species.

        Parameters
        ----------
        wdot      : numpy array, shape (44,), modified in-place.
        sc        : numpy array, shape (44,), molar concentrations [mol/m³].
        T         : float, gas temperature [K].
        Te        : float, electron temperature [K].
        EN        : float, reduced electric field (unused stub).
        ener_exch : numpy array, shape (44,), energy exchange, modified in-place.
        """
        # --- temperature factors ---
        invT = 1.0 / T
        logT = math.log(T / 300.0)   # log(T/300), same as tc[0] in C++

        # reference concentration: P_atm / (R T)  [mol/m³]
        refC    = 101325.0 / 8.31446 * invT
        refCinv = 1.0 / refC

        # --- initialise wdot ---
        for i in range(NUM_SPECIES):
            wdot[i] = 0.0

        # --- mixture concentration ---
        mixture = float(np.sum(sc))

        # --- Gibbs free energies ---
        g_RT = np.zeros(NUM_SPECIES)
        gibbs(g_RT, T)

        # --- electron energy terms ---
        ne = sc[0] * NA          # electron number density [1/m³]
        Ue = 1.5 * Te * ne * 1.380649e-23

        # --- electron temperature quantities ---
        if Te <= 0.0:
            invTe   = 1.0
            TeeV    = 1e-30
            logTe   = math.log(1e-30)
            invTeeV = 1.0 / 1e-30
        else:
            invTe   = 1.0 / Te
            TeeV    = Te / 11595.0
            logTe   = math.log(TeeV)
            invTeeV = 1.0 / TeeV

        Te_pow    = [logTe ** j for j in range(9)]
        invTe_pow = [invTeeV ** (j + 1) for j in range(4)]

        # --- iterate over reactions ---
        for rxn in self.reactions:
            self._apply_reaction(
                rxn, wdot, sc, T, Te, TeeV, logT, invT,
                logTe, invTeeV, Te_pow, invTe_pow,
                mixture, g_RT, refC, refCinv, ener_exch, Ue
            )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _compute_k_f(self, rxn, TeeV, logT, invT,
                     logTe, Te_pow, invTe_pow, invTeeV):
        """Evaluate the forward rate constant for a single reaction."""
        rt = rxn['rate_type']

        if rt == 'Jfit':
            coefs = rxn['Jfit_coefs']
            A     = rxn['Jfit_A']
            s     = sum(coefs[j] * Te_pow[j] for j in range(len(coefs)))
            k_f   = A * math.exp(s) * NA
            k_max = rxn.get('k_max')
            if k_max is not None:
                k_f = min(k_f, k_max)
            return k_f

        if rt == 'Ffit':
            coefs = rxn['Ffit_coefs']
            A     = rxn['Ffit_A']
            s     = sum(coefs[j] * invTe_pow[j] for j in range(len(coefs)))
            return A * math.exp(s) * NA

        if rt == 'conditional':
            threshold = rxn['TeeV_threshold']
            if TeeV < threshold:
                coefs = rxn['low_Jfit_coefs']
                A     = rxn['low_Jfit_A']
                s     = sum(coefs[j] * Te_pow[j] for j in range(len(coefs)))
                return A * math.exp(s) * NA
            else:
                coefs = rxn['high_Ffit_coefs']
                A     = rxn['high_Ffit_A']
                s     = sum(coefs[j] * invTe_pow[j] for j in range(len(coefs)))
                return A * math.exp(s) * NA

        if rt == 'Arrhenius':
            A    = rxn['A']
            beta = rxn.get('beta', 0.0)
            Ea   = rxn.get('Ea', 0.0)
            return A * math.exp(beta * logT - Ea * invT)

        if rt == 'constant':
            return rxn['k_f']

        raise ValueError(f"Unknown rate_type '{rt}' in reaction id={rxn.get('id')}")

    def _compute_qr(self, rxn, k_f, g_RT, refC, refCinv, mixture):
        """
        Compute the reverse rate contribution.

        NOTE: In the current prodfunc.cpp all reverse rates are multiplied by
        (0.0), so qr is always zero.  This function faithfully replicates that
        behaviour while still computing the Gibbs-based prefactor for
        completeness.
        """
        pairs = rxn.get('gibbs_exponent', [])
        if not pairs:
            return 0.0

        gibbs_sum = sum(coef * g_RT[idx] for idx, coef in pairs)
        rf = rxn.get('refC_factor', 0)
        if rf == 1:
            refC_val = refC
        elif rf == -1:
            refC_val = refCinv
        else:
            refC_val = 1.0

        Corr = mixture if rxn.get('has_third_body') else 1.0
        # Always zero because the C++ code multiplies by (0.0)
        return Corr * k_f * math.exp(-gibbs_sum) * refC_val * 0.0

    def _apply_reaction(self, rxn, wdot, sc, T, Te, TeeV, logT, invT,
                        logTe, invTeeV, Te_pow, invTe_pow,
                        mixture, g_RT, refC, refCinv, ener_exch, Ue):
        """Apply a single reaction to wdot and ener_exch."""
        k_f = self._compute_k_f(rxn, TeeV, logT, invT,
                                 logTe, Te_pow, invTe_pow, invTeeV)

        # Forward concentration product
        qf_conc = 1.0
        for idx in rxn.get('reactants', []):
            qf_conc *= sc[idx]
        if rxn.get('has_third_body'):
            qf_conc *= mixture

        qf = k_f * qf_conc
        qr = self._compute_qr(rxn, k_f, g_RT, refC, refCinv, mixture)

        qdot = qf - qr

        # Update species production rates
        for sp_idx, coeff in rxn.get('stoichiometry', []):
            wdot[sp_idx] += coeff * qdot

        # Energy exchange
        ee = rxn.get('energy_exchange')
        if ee is not None:
            comp_ener_exch(
                qf, qr, sc, k_f,
                ee['rxntype'], ee['eexci'], ee['elidx'],
                ener_exch, Ue, T, Te
            )


# ---------------------------------------------------------------------------
# CLI smoke test
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    reactions_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reactions')
    db = ReactionDB(reactions_dir)
    print(f"Loaded {len(db.reactions)} reactions")

    # Quick sanity check with dummy inputs
    wdot      = np.zeros(NUM_SPECIES)
    sc        = np.zeros(NUM_SPECIES)
    ener_exch = np.zeros(NUM_SPECIES)

    # Set a small amount of each species so no division-by-zero
    sc[:] = 1e-6

    try:
        db.production_rate(wdot, sc, T=1000.0, Te=10000.0, EN=0.0,
                           ener_exch=ener_exch)
        n_nonzero = int(np.count_nonzero(wdot))
        print(f"production_rate() ran successfully; {n_nonzero} non-zero wdot entries")
    except Exception as exc:
        import traceback
        traceback.print_exc()
        raise
