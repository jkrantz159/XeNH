"""
ParallelXeModel - Xenon isotope evolution model for Earth's mantle

Models the evolution of 128Xe and 130Xe in Earth's mantle over
geological time using a box model with sigmoidal downwelling growth.

Based on: Parai and Mukhopadhyay (2018)
"""

import numpy as np
from model_config import ModelConfig
from model_utils import ModelUtils


def parallel_xe_model(xe_cap, alpha, beta, eta, resfrac, lv_frac,
                     plot_check, atm, t, T, count=0):
    """
    Xenon isotope evolution model

    Parameters:
        xe_cap (float): Xe carrying capacity (atoms/gram)
        alpha (float): Growth rate parameter (/yr)
        beta (float): Sigmoid inflection point (years)
        eta (float): Processing rate parameter (/yr)
        resfrac (float): Reservoir fraction
        lv_frac (float): Late veneer fraction (%)
        plot_check (int): 1 to generate plots, 0 to skip
        atm (ndarray): Atmospheric 128Xe/130Xe evolution
        t (ndarray): Time vector (years)
        T (float): Total age of Earth (years)
        count (int): Iteration counter

    Returns:
        int: 1 if success, 0 if failure
    """
    # Initialize arrays
    xe130_m = np.zeros(len(t))
    xe128_m = np.zeros(len(t))

    # Physical constants
    mres = ModelConfig.EARTH_MASS_GRAMS * resfrac
    me = ModelConfig.EARTH_MASS_GRAMS

    # Calculate initial mantle composition
    mol130_m = (ModelConfig.XE130_CONCENTRATION_CC * me) * (lv_frac / 100)
    initial130 = mol130_m * ModelConfig.AVOGADRO_NUMBER / mres
    initial128 = initial130 * ModelConfig.XE128_130_RATIO_AVCC

    # Calculate sigmoidal downwelling evolution
    xed = ModelUtils.calculate_sigmoidal_downwelling(xe_cap, alpha, beta, t)

    # Time evolution loop
    for i in range(1, len(t)):
        # Calculate mass processed
        dM = ModelUtils.calculate_mass_processed(
            eta, t[i], t[i-1], T, ModelConfig.PRESENT_DAY_PROCESSING_RATE)
        dM_mres = dM / mres

        # Get previous values
        if i == 1:
            xe130_mlast = initial130
            xe128_mlast = initial128
            atmlast = atm[0]
            xedlast = xed[i]
        else:
            xe130_mlast = xe130_m[i-1]
            xe128_mlast = xe128_m[i-1]
            atmlast = atm[i-1]
            xedlast = xed[i-1]

        # Update concentrations
        xe130_m[i] = ModelUtils.update_concentration(xe130_mlast, dM_mres, xedlast, 1)
        xe128_m[i] = ModelUtils.update_concentration(xe128_mlast, dM_mres, xedlast, atmlast)

    # Check success criteria
    succ1 = ModelUtils.check_success_criteria(
        xe130_m[-1], ModelConfig.XE130_MANTLE_MIN, ModelConfig.XE130_MANTLE_MAX)

    ratio = xe128_m[-1] / xe130_m[-1]
    succ2 = ModelUtils.check_success_criteria(
        ratio, ModelConfig.XE128_130_MANTLE_MIN, ModelConfig.XE128_130_MANTLE_MAX)

    return 1 if (succ1 and succ2) else 0
