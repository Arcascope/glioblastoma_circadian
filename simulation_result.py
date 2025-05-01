import numpy as np


class SimulationResult:
    def __init__(self,
                 sol,
                 dlmos,
                 zeitgeber,
                 pills_local_time,
                 pills_corrected_for_dst,
                 pill_melatonin_offsets,
                 actogram):
        self.sol = sol
        self.dlmos = dlmos
        self.zeitgeber = zeitgeber
        self.pills_local_time = np.array(pills_local_time)
        self.pills_corrected_for_dst = np.array(pills_corrected_for_dst)
        self.pill_melatonin_offsets = np.array(pill_melatonin_offsets)
        self.actogram = actogram
