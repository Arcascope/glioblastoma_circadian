import pickle

import matplotlib.pyplot as plt
import utils
import numpy as np


def make_histograms(participant_array, circular=False):
    zeitgeber_times_am = []
    zeitgeber_times_pm = []
    circadian_times_am = []
    circadian_times_pm = []

    for subject in participant_array:
        if subject.group == utils.Group.AM:
            zeitgeber_times_am += subject.all_pills_local_time.tolist()
            circadian_times_am += subject.all_pill_offsets.tolist()
        else:
            zeitgeber_times_pm += subject.all_pills_local_time.tolist()
            circadian_times_pm += subject.all_pill_offsets.tolist()

    zeitgeber_times_am = np.array(zeitgeber_times_am)
    zeitgeber_times_pm = np.array(zeitgeber_times_pm)

    circadian_times_am = np.array(circadian_times_am)
    circadian_times_pm = np.array(circadian_times_pm)

    if circular:
        fig, axes = plt.subplots(1, 2, figsize=(
            12, 4), subplot_kw={'projection': 'polar'})
        num_bins = 24  # Number of bins (intervals)

        alpha = 1.0
        axes[0].hist(np.radians(zeitgeber_times_am * 15), bins=num_bins, color=utils.AM_COLOR, alpha=alpha,
                     label='Wall Clock AM')
        axes[0].hist(np.radians(zeitgeber_times_pm * 15), bins=num_bins, color=utils.PM_COLOR, alpha=alpha,
                     label='Wall Clock PM')
        axes[1].hist(np.radians(circadian_times_am * 15), bins=num_bins, color=utils.AM_COLOR, alpha=alpha,
                     label='Biological Time AM')
        axes[1].hist(np.radians(circadian_times_pm * 15), bins=num_bins, color=utils.PM_COLOR, alpha=alpha,
                     label='Biological Time PM')

        tick_interval = 360 // 24
        num_ticks = 24

        x_tick_labels = [str(i) for i in range(0, 24, 24 // num_ticks)]

        for index in range(2):
            axes[index].set_xticks(np.radians(
                range(0, 360, tick_interval)))  # Set tick locations
            axes[index].set_xticklabels(x_tick_labels)  # Set tick labels
            # Set zero angle at the top (North)
            axes[index].set_theta_zero_location('N')
            axes[index].set_yticklabels([])
            axes[index].legend(loc='upper right')
            axes[index].xaxis.grid(linewidth=0.2, color='gray', alpha=0.5)
            axes[index].yaxis.grid(linewidth=0.2, color='gray', alpha=0.5)

        axes[0].set_title('Wall clock dosing (AM vs PM)')
        axes[1].set_title('Biological time dosing (AM vs PM)')
    else:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        alpha = 0.9  # Adjust alpha for transparency
        zt_cut_point = 15
        ct_cut_point = 8

        zeitgeber_times_am[zeitgeber_times_am > zt_cut_point] -= 24
        zeitgeber_times_pm[zeitgeber_times_pm > zt_cut_point] -= 24

        circadian_times_am[circadian_times_am > ct_cut_point] -= 24
        circadian_times_pm[circadian_times_pm > ct_cut_point] -= 24

        # Change to "hours after DLMO"
        circadian_times_am = - circadian_times_am
        circadian_times_pm = - circadian_times_pm

        bio_start = -6
        bio_end = bio_start + 24

        axes[0].hist(zeitgeber_times_am, bins=np.arange(-8, 16, 0.5),
                     color=utils.AM_COLOR, alpha=alpha, label='Wall Clock AM')
        axes[0].hist(zeitgeber_times_pm, bins=np.arange(-8, 16, 0.5),
                     color=utils.PM_COLOR, alpha=alpha, label='Wall Clock PM')

        axes[1].hist(circadian_times_am, bins=np.arange(bio_start, bio_end, 0.5),
                     color=utils.AM_COLOR, alpha=alpha, label='Biological Time AM')
        axes[1].hist(circadian_times_pm, bins=np.arange(bio_start, bio_end, 0.5),
                     color=utils.PM_COLOR, alpha=alpha, label='Biological Time PM')

        max_circadian_pm = np.max(circadian_times_pm)
        within_4_hours_am = np.sum(
            np.abs(circadian_times_am - max_circadian_pm) <= 4)
        percent_within_4_hours_am = (
            within_4_hours_am / len(circadian_times_am)) * 100
        print(
            f"Number of circadian_times_am within 4 hours of the highest circadian_times_pm value: {within_4_hours_am} ({percent_within_4_hours_am:.2f}%)")

        min_circadian_am = np.min(circadian_times_am)
        within_4_hours_pm = np.sum(
            np.abs(circadian_times_pm - min_circadian_am) <= 4)
        percent_within_4_hours_pm = (
            within_4_hours_pm / len(circadian_times_pm)) * 100
        print(
            f"Number of circadian_times_pm within 4 hours of the lowest circadian_times_am value: {within_4_hours_pm} ({percent_within_4_hours_pm:.2f}%)")

        max_zeitgeber_pm = np.max(zeitgeber_times_pm)
        within_4_hours_am_zt = np.sum(
            np.abs(zeitgeber_times_am - max_zeitgeber_pm) <= 4)
        percent_within_4_hours_am_zt = (
            within_4_hours_am_zt / len(zeitgeber_times_am)) * 100
        print(
            f"Number of zeitgeber_times_am within 4 hours of the highest zeitgeber_times_pm value: {within_4_hours_am_zt} ({percent_within_4_hours_am_zt:.2f}%)")

        min_zeitgeber_am = np.min(zeitgeber_times_am)
        within_4_hours_pm_zt = np.sum(
            np.abs(zeitgeber_times_pm - min_zeitgeber_am) <= 4)
        percent_within_4_hours_pm_zt = (
            within_4_hours_pm_zt / len(zeitgeber_times_pm)) * 100
        print(
            f"Number of zeitgeber_times_pm within 4 hours of the lowest zeitgeber_times_am value: {within_4_hours_pm_zt} ({percent_within_4_hours_pm_zt:.2f}%)")
        spread_am = np.max(circadian_times_am) - np.min(circadian_times_am)
        spread_pm = np.max(circadian_times_pm) - np.min(circadian_times_pm)
        print(f"Spread of circadian_times_am: {spread_am}")
        print(f"Spread of circadian_times_pm: {spread_pm}")

        spread_zeitgeber_am = np.max(
            zeitgeber_times_am) - np.min(zeitgeber_times_am)
        spread_zeitgeber_pm = np.max(
            zeitgeber_times_pm) - np.min(zeitgeber_times_pm)
        print(f"Spread of zeitgeber_times_am: {spread_zeitgeber_am}")
        print(f"Spread of zeitgeber_times_pm: {spread_zeitgeber_pm}")

        plt.xticks(np.arange(-8, 16, 4))
        current_xticks = axes[0].get_xticks()
        hour_labels = ["{:02d}:00".format(round(np.mod(hour + 24, 24))) for hour in
                       current_xticks]
        label_size = 14
        axes[0].set_xticks(current_xticks, hour_labels)
        axes[0].set_xlim([-8, 16])
        axes[0].set_xlabel('Wall Clock Time', fontsize=label_size)
        axes[0].legend(loc='upper left')
        axes[0].set_ylabel('Count', fontsize=label_size)

        axes[1].set_xticks(range(bio_start, bio_end, 2))
        axes[1].legend(loc='upper left')

        current_xticks_bio = axes[1].get_xticks()
        bio_hour_labels = [np.mod(hour, 24) for hour in current_xticks_bio]
        axes[1].set_xticks(current_xticks_bio, bio_hour_labels)
        axes[1].set_xlabel(
            'Biological Time (hours after DLMO mod 24)', fontsize=label_size)

        axes[1].set_xlim([bio_start, bio_end])

        for ax in axes.flatten():
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)

        plt.tight_layout()

    plt.savefig("figures/histogram.png", dpi=utils.DPI)
    plt.show()


if __name__ == '__main__':
    with open(utils.snapshot_path, 'rb') as input_file:
        people = pickle.load(input_file)
    make_histograms(people)
