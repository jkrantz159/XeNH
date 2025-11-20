"""
ModelConfig - Configuration and constants for XeNH mantle recycling models

This module contains all physical constants, initial conditions, and
model parameters used in the H-N-Xe mantle recycling simulations.

References:
- Marty (2012): Noble gas concentrations in chondrites
- Pepin (2000): Isotopic composition of primordial Xe
- Parai and Mukhopadhyay (2018): Model framework
- Barry and Hilton (2016): Nitrogen isotope systematics
- Porcelli, Ballentine, and Wieler (2002): Earth atmosphere composition
"""

import numpy as np


class ModelConfig:
    """Configuration and constants for XeNH models"""

    # Physical Constants
    EARTH_MASS_GRAMS = 5.972e27  # Mass of Earth (g)
    EARTH_AGE_YEARS = 4.568e9    # Age of Earth (years)
    AVOGADRO_NUMBER = 6.022e23   # Avogadro's number (atoms/mol)

    # Mantle Properties
    CONVECTING_MANTLE_FRACTION = 0.9       # 90% of mantle is convecting
    CONVECTING_MANTLE_MASS = 3.6e27        # Mass of convecting mantle (g)

    # Mantle Processing Rates
    PRESENT_DAY_PROCESSING_RATE = 6.1e17   # g/yr, based on He flux at ridges

    # Late Veneer
    LATE_VENEER_FRACTION_DEFAULT = 1.0     # 1% of Earth mass

    # Xenon Initial Conditions (AVCC)
    XE130_CONCENTRATION_CC = 5.38e-14      # mol/g (Marty 2012)
    XE128_130_RATIO_AVCC = 0.5073          # Pepin 2000

    # Xenon Present-Day Atmosphere
    XE128_132_RATIO_ATM = 0.0714
    XE130_132_RATIO_ATM = 0.1514
    XE128_130_RATIO_ATM_TODAY = 0.4716     # 0.0714/0.1514
    XE128_130_RATIO_ATM_INITIAL = 0.5178   # -39 per mille fractionation
    ATMOSPHERE_EVOLUTION_TIME = 2.568e9    # years

    # Xenon Success Criteria
    XE130_MANTLE_MIN = 4.3e5
    XE130_MANTLE_MAX = 9.2e5
    XE128_130_MANTLE_MIN = 0.475
    XE128_130_MANTLE_MAX = 0.478

    # Nitrogen Initial Conditions
    N_CONCENTRATION_CC = 0.001519          # weight fraction (Sephton et al. 2003)
    N14_FRACTION = 0.99636
    N15_14_RATIO_INITIAL = 2.3e-3          # Owen et al. 2001
    N15_14_RATIO_SEDIMENT = 0.003671565    # delta-15N = +5 per mille
    N_ATOMIC_MASS_14 = 14

    # Nitrogen Success Criteria
    N14_MANTLE_MIN_MOL = 7.06e19
    N14_MANTLE_MAX_MOL = 9.78e21
    N15_14_RATIO_MANTLE_MIN = 0.0036275
    N15_14_RATIO_MANTLE_MAX = 0.0036425

    # Neon Initial Conditions
    NE22_CONCENTRATION_CC = 1.62e-12       # mol/g (Marty 2012)
    NE20_22_RATIO_AVCC_PSN = 13.36         # Williams and Mukhopadhyay 2018

    # Neon Success Criteria
    NE22_MANTLE_CENTER = 5.8e-15
    NE22_MANTLE_UNCERTAINTY = 3.2e-15
    NE20_22_RATIO_MANTLE_CENTER = 12
    NE20_22_RATIO_MANTLE_TOLERANCE = 0.5
    NE130_FRACTIONATION = 9.80

    # Model Parameters - Default Ranges
    ALPHA_MIN = 1e-10
    ALPHA_MAX = 1e-7
    BETA_MIN = 0
    BETA_MAX = 10e9
    ETA_MIN = 7.5e-10
    ETA_MAX = 8.0e-10
    XE_CAPACITY_MIN = 5
    XE_CAPACITY_MAX = 5e8
    N_CAPACITY_MIN = 4
    N_CAPACITY_MAX = 4e17

    # Time Vector Parameters
    TIME_START = 0
    TIME_END = 4.568e9
    TIME_STEP_EARLY = 0.1e6
    TIME_STEP_MIDDLE = 1e6
    TIME_STEP_LATE = 5e6
    TIME_EARLY_END = 200e6
    TIME_MIDDLE_START = 201e6
    TIME_MIDDLE_END = 3300e6
    TIME_LATE_START = 3305e6

    @staticmethod
    def create_time_vector():
        """Create the time vector for model integration"""
        t = np.concatenate([
            np.arange(ModelConfig.TIME_START,
                     ModelConfig.TIME_EARLY_END + ModelConfig.TIME_STEP_EARLY,
                     ModelConfig.TIME_STEP_EARLY),
            np.arange(ModelConfig.TIME_MIDDLE_START,
                     ModelConfig.TIME_MIDDLE_END + ModelConfig.TIME_STEP_MIDDLE,
                     ModelConfig.TIME_STEP_MIDDLE),
            np.arange(ModelConfig.TIME_LATE_START,
                     ModelConfig.TIME_END + ModelConfig.TIME_STEP_LATE,
                     ModelConfig.TIME_STEP_LATE),
            [ModelConfig.TIME_END + 3e6]
        ])
        return t

    @staticmethod
    def create_atmosphere_evolution(t):
        """Create 128Xe/130Xe atmospheric evolution"""
        atm = np.zeros(len(t))
        for i, time in enumerate(t):
            if time < ModelConfig.ATMOSPHERE_EVOLUTION_TIME:
                # Linear evolution from initial to modern
                atm[i] = (ModelConfig.XE128_130_RATIO_ATM_INITIAL -
                         (time / ModelConfig.ATMOSPHERE_EVOLUTION_TIME) *
                         (ModelConfig.XE128_130_RATIO_ATM_INITIAL -
                          ModelConfig.XE128_130_RATIO_ATM_TODAY))
            else:
                atm[i] = ModelConfig.XE128_130_RATIO_ATM_TODAY
        return atm

    @staticmethod
    def ratio_to_delta_xe(ratio):
        """Convert 128Xe/130Xe ratio to delta notation (per mille)"""
        return ((ratio / ModelConfig.XE128_130_RATIO_ATM_TODAY) - 1) * 1000

    @staticmethod
    def ratio_to_delta_n(ratio):
        """Convert 15N/14N ratio to delta notation (per mille)"""
        air_ratio = (1 - ModelConfig.N14_FRACTION) / ModelConfig.N14_FRACTION
        return ((ratio / air_ratio) - 1) * 1000
