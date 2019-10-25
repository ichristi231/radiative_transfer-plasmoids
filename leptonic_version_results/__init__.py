from pathlib import Path
import numpy as np

tables_dir = Path.cwd() / "leptonic_version_results"

electron_distribution = np.loadtxt(
    tables_dir / "electron_distribution_plasmoid_1.txt")

co_moving_lumionsity = np.loadtxt(
    tables_dir / "co_moving_nu_l_nu_plasmoid_1.txt")

electron_Lorentz_factor = np.loadtxt(
    tables_dir / "log_10_electron_Lorentz_factor_array.txt")

photon_frequency = np.loadtxt(
    tables_dir / "log_10_photon_frequency_array.txt")
