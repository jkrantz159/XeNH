"""
ParallelNewNModel - Nitrogen isotope evolution model for Earth's mantle

Models the evolution of 14N and 15N in Earth's mantle over
geological time using a box model with sigmoidal downwelling growth.

Based on: Parai and Mukhopadhyay (2018) and Barry and Hilton (2016)
"""

import numpy as np
from model_config import ModelConfig
from model_utils import ModelUtils


def parallel_n_model(n_cap, alpha, beta, eta, resfrac, lv_frac,
                    plot_check, t, T, count=0):
    """
    Nitrogen isotope evolution model

    Parameters:
        n_cap (float): N carrying capacity (atoms/gram)
        alpha (float): Growth rate parameter (/yr)
        beta (float): Sigmoid inflection point (years)
        eta (float): Processing rate parameter (/yr)
        resfrac (float): Reservoir fraction
        lv_frac (float): Late veneer fraction (%)
        plot_check (int): 1 to generate plots, 0 to skip
        t (ndarray): Time vector (years)
        T (float): Total age of Earth (years)
        count (int): Iteration counter

    Returns:
        int: 1 if success, 0 if failure
    """
    # Initialize arrays
    n14_m = np.zeros(len(t))
    n15_m = np.zeros(len(t))

    # Physical constants
    mres = ModelConfig.EARTH_MASS_GRAMS * resfrac
    me = ModelConfig.EARTH_MASS_GRAMS

    # Calculate initial mantle composition
    g14_m = (ModelConfig.N_CONCENTRATION_CC * me * lv_frac / 100)
    initial14 = (g14_m / ModelConfig.N_ATOMIC_MASS_14 * ModelConfig.N14_FRACTION *
                ModelConfig.AVOGADRO_NUMBER) / mres
    initial15 = initial14 * ModelConfig.N15_14_RATIO_INITIAL

    # Calculate sigmoidal downwelling evolution
    nd = ModelUtils.calculate_sigmoidal_downwelling(n_cap, alpha, beta, t)

    # Time evolution loop
    for i in range(1, len(t)):
        # Calculate mass processed
        dM = ModelUtils.calculate_mass_processed(
            eta, t[i], t[i-1], T, ModelConfig.PRESENT_DAY_PROCESSING_RATE)
        dM_mres = dM / mres

        # Get previous values
        if i == 1:
            n14_mlast = initial14
            n15_mlast = initial15
            ndlast = nd[i]
        else:
            n14_mlast = n14_m[i-1]
            n15_mlast = n15_m[i-1]
            ndlast = nd[i-1]

        # Update concentrations
        n14_m[i] = ModelUtils.update_concentration(n14_mlast, dM_mres, ndlast, 1)
        n15_m[i] = ModelUtils.update_concentration(
            n15_mlast, dM_mres, ndlast, ModelConfig.N15_14_RATIO_SEDIMENT)

    # Check success criteria
    n14_min = (ModelConfig.N14_MANTLE_MIN_MOL * ModelConfig.N14_FRACTION *
              ModelConfig.AVOGADRO_NUMBER / mres)
    n14_max = (ModelConfig.N14_MANTLE_MAX_MOL * ModelConfig.N14_FRACTION *
              ModelConfig.AVOGADRO_NUMBER / mres)
    succ1 = ModelUtils.check_success_criteria(n14_m[-1], n14_min, n14_max)

    ratio = n15_m[-1] / n14_m[-1]
    succ2 = ModelUtils.check_success_criteria(
        ratio, ModelConfig.N15_14_RATIO_MANTLE_MIN, ModelConfig.N15_14_RATIO_MANTLE_MAX)

    return 1 if (succ1 and succ2) else 0
