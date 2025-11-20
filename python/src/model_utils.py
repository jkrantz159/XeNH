"""
ModelUtils - Utility functions for XeNH mantle recycling models

This module contains common calculations used across different
geochemical models to reduce code duplication.
"""

import numpy as np
import os
from pathlib import Path


class ModelUtils:
    """Utility functions for geochemical models"""

    @staticmethod
    def calculate_mass_processed(eta, t_curr, t_prev, T, Qp):
        """
        Calculate mass of mantle processed in time step

        Parameters:
            eta (float): Processing rate parameter (/yr)
            t_curr (float): Current time (years)
            t_prev (float): Previous time (years)
            T (float): Total age of Earth (years)
            Qp (float): Present day processing rate (g/yr)

        Returns:
            float: Mass of mantle processed (g)
        """
        dM = (Qp / eta) * (np.exp(eta * (T - t_prev)) - np.exp(eta * (T - t_curr)))
        return dM

    @staticmethod
    def calculate_sigmoidal_downwelling(capacity, alpha, beta, t):
        """
        Calculate downwelling using sigmoidal model

        Parameters:
            capacity (float): Carrying capacity
            alpha (float): Growth rate parameter (/yr)
            beta (float): Inflection point (years)
            t (ndarray): Time vector (years)

        Returns:
            ndarray: Downwelling flux at each time
        """
        downwelling = capacity / (1 + np.exp(-alpha * (t - beta)))
        return downwelling

    @staticmethod
    def update_concentration(prev_conc, dM_Mres, downwelling, downwelling_ratio=1.0):
        """
        Update isotope concentration using box model

        Parameters:
            prev_conc (float): Previous concentration in mantle
            dM_Mres (float): Normalized mass processed (dM/Mres)
            downwelling (float): Downwelling concentration
            downwelling_ratio (float): Isotopic ratio of downwelling

        Returns:
            float: Updated concentration
        """
        concentration = prev_conc + dM_Mres * (downwelling * downwelling_ratio - prev_conc)
        return concentration

    @staticmethod
    def check_success_criteria(value, min_val, max_val):
        """Check if value falls within acceptable range"""
        return (value >= min_val) and (value <= max_val)

    @staticmethod
    def ensure_directory_exists(dir_path):
        """Create directory if it doesn't exist"""
        Path(dir_path).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def get_output_path(filename, subdir='results'):
        """Construct output file path"""
        project_root = Path(__file__).parent.parent
        output_dir = project_root / subdir
        output_dir.mkdir(parents=True, exist_ok=True)
        return output_dir / filename
