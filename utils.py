from datetime import datetime

import numpy as np
import pandas as pd
from enum import Enum
import matplotlib

# Set the font globally to Helvetica
matplotlib.rcParams['font.family'] = 'Helvetica'

output_folder = "debug_outputs/"
MODE = "activity"  # "light+activity"  # activity, light <-- Three choices for model zeitgeber, story is the same for all of them

snapshot_path = f"snapshots/{MODE}_saved_people.pkl"
DPI = 300
LIGHT_SCALAR = 1.0 if MODE == "light" else 0.2

TEAL_REGULAR = (54 / 255, 173 / 255, 166 / 255)
TEAL_LIGHT = (197 / 255, 238 / 255, 236 / 255)
ORANGE_REGULAR = (252 / 255, 136 / 255, 64 / 255)
ORANGE_LIGHT = (255 / 255, 226 / 255, 208 / 255)
GOLD_REGULAR = (240 / 255, 196 / 255, 74 / 255)
GOLD_LIGHT = (255 / 255, 239 / 255, 195 / 255)
YELLOW_REGULAR = (239 / 255, 221 / 255, 24 / 255)
YELLOW_LIGHT = (255 / 255, 253 / 255, 227 / 255)
GREEN_REGULAR = (137 / 255, 166 / 255, 100 / 255)
GREEN_LIGHT = (216 / 255, 230 / 255, 197 / 255)
PERIWINKLE_REGULAR = (133 / 255, 148 / 255, 214 / 255)
PERIWINKLE_LIGHT = (217 / 255, 223 / 255, 255 / 255)
PURPLE_REGULAR = (129 / 255, 59 / 255, 140 / 255)
PURPLE_LIGHT = (252 / 255, 230 / 255, 255 / 255)
NEUTRAL_WHITE = (255 / 255, 255 / 255, 255 / 255)
NEUTRAL_LINK_BLUE = (43 / 255, 64 / 255, 144 / 255)
NEUTRAL_LINK_BLUE_LIGHT = (153 / 255, 200 / 255, 255 / 255)
NEUTRAL_GRAY_5 = (48 / 255, 47 / 255, 46 / 255)
NEUTRAL_GRAY_4 = (71 / 255, 69 / 255, 64 / 255)
NEUTRAL_GRAY_3 = (94 / 255, 92 / 255, 92 / 255)
NEUTRAL_GRAY_2 = (237 / 255, 237 / 255, 237 / 255)
NEUTRAL_GRAY_1 = (248 / 255, 248 / 255, 248 / 255)
RED_REGULAR = (224 / 255, 80 / 255, 67 / 255)
RED_LIGHT = (251 / 255, 188 / 255, 183 / 255)

AM_COLOR = YELLOW_REGULAR
PM_COLOR = TEAL_REGULAR
BACKGROUND_COLOR = NEUTRAL_GRAY_2

# Map for subject IDs to numbers for paper
subject_id_to_number = {
    "P33": 1,
    "P34": 2,
    "P35": 3,
    "P36": 4,
    "P38": 5,
    "P39": 6,
    "P41": 7,
    "P42": 8,
    "P43": 9,
    "P44": 10
}

END_OF_TRIAL = datetime(day=1, month=8, year=2022, hour=0,
                        minute=0, second=0, microsecond=0)


def circular_mean(radians):
    sin_sum = np.sum(np.sin(radians))
    cos_sum = np.sum(np.cos(radians))
    mean_angle = np.arctan2(sin_sum, cos_sum)
    return mean_angle


def read_outcomes_spreadsheet():
    file_path = 'data/Chronotherapy Results - Updated based on conversation with meeting with Omar 10-3.xlsx'
    data = pd.read_excel(file_path)

    # Preprocess the data
    data['Date of Death'] = pd.to_datetime(
        data['Date of Death']).fillna(pd.Timestamp(END_OF_TRIAL))
    data['Surgery date'] = pd.to_datetime(data['Surgery date'])
    data['Start Date of Drug'] = pd.to_datetime(data['Start Date of Drug'])
    data['End Date of Drug'] = pd.to_datetime(data['End Date of Drug'])
    data['Date of Pseudoprogression'] = pd.to_datetime(
        data['Date of Pseudoprogression'])
    data['Date of Progression (1)'] = pd.to_datetime(
        data['Date of Progression (1)'])
    data['Date of Progression (2)'] = pd.to_datetime(
        data['Date of Progression (2)'])
    data['Alive'] = data['Date of Death'] == pd.Timestamp(END_OF_TRIAL)

    return data


def phase_ic_guess(time_of_day: float):
    time_of_day = np.fmod(time_of_day, 24.0)

    # Wake at 7 am after 8 hours of sleep, state at 00:00
    psi = 1.65238233

    # Convert to radians, add to phase
    psi += time_of_day * np.pi / 12
    return psi


class Group(Enum):
    AM = "AM"
    PM = "PM"

    def color(self):
        if self.value == "AM":
            return AM_COLOR
        elif self.value == "PM":
            return PM_COLOR
