from simulation_result import SimulationResult
import numpy as np
import pandas as pd
from scipy.stats import circstd

from utils import Group


class PersonMetrics:
    def __init__(self, subject_id: str, group: Group, simulation_result: SimulationResult):
        self.subject_id = subject_id
        self.group = group
        self.shifted_melatonin_onsets = np.mod(np.array(simulation_result.dlmos) - 12, 24)
        self.mean_melatonin_onset = np.nanmean(self.shifted_melatonin_onsets) + 12
        self.std_melatonin_onset = np.nanstd(self.shifted_melatonin_onsets)
        self.moving_std = pd.Series(self.shifted_melatonin_onsets).rolling(7).std().values
        self.moving_std_mean = np.nanmean(self.moving_std)
        self.pill_std = circstd(simulation_result.pills_local_time / 24 * 2 * np.pi, nan_policy='omit') * 24 / (2 * np.pi)
        self.pill_offset_std = circstd(simulation_result.pill_melatonin_offsets / 24 * 2 * np.pi, nan_policy='omit') * 24 / (
                2 * np.pi)
        self.all_pill_offsets = simulation_result.pill_melatonin_offsets
        self.all_pills_local_time = simulation_result.pills_local_time
        self.actogram = simulation_result.actogram

    def moving_std_with_window(self, window):
        values_to_calculate = np.array(self.shifted_melatonin_onsets) / 24 * 2 * np.pi
        return pd.Series(values_to_calculate).rolling(window).apply(circstd).values * 24 / (2 * np.pi)

    def average_of_tail_with_length(self, number_to_use, window_size=3):
        return np.mean(self.moving_std_with_window(window_size)[-number_to_use:])
