import pickle
import numpy as np
import matplotlib.pyplot as plt
import utils


def make_compliance_figure(participant_array):
    zeitgeber_times_am = []
    zeitgeber_times_pm = []
    circadian_times_am = []
    circadian_times_pm = []

    all_zt_stds = []
    all_ct_stds = []

    for subject in participant_array:
        if subject.group == utils.Group.AM:
            zeitgeber_times_am.append(subject.pill_std)
            circadian_times_am.append(subject.pill_offset_std)
        else:
            zeitgeber_times_pm.append(subject.pill_std)
            circadian_times_pm.append(subject.pill_offset_std)

        all_zt_stds.append(subject.pill_std)
        all_ct_stds.append(subject.pill_offset_std)

    zeitgeber_times_am = np.array(zeitgeber_times_am)
    zeitgeber_times_pm = np.array(zeitgeber_times_pm)
    circadian_times_am = np.array(circadian_times_am)
    circadian_times_pm = np.array(circadian_times_pm)
    all_zt_stds = np.array(all_zt_stds)
    all_ct_stds = np.array(all_ct_stds)

    fig, axes = plt.subplots(1, 1, figsize=(12, 4))

    y_bar_height = 1.8
    axes.plot([0, 1], [y_bar_height, y_bar_height], color='black')
    axes.text(0.5, y_bar_height * 0.95, '*', ha='center',
              va='bottom', color='black', fontsize=44)

    axes.plot([2, 3], [y_bar_height, y_bar_height], color='black')
    axes.text(2.5,  y_bar_height * 0.95, '*', ha='center',
              va='bottom', color='black', fontsize=44)
    axes.bar([0, 1, 2, 3], [np.nanmean(zeitgeber_times_am),
                            np.nanmean(circadian_times_am),
                            np.nanmean(zeitgeber_times_pm),
                            np.nanmean(circadian_times_pm)
                            ],
             yerr=[np.nanstd(zeitgeber_times_am),
                   np.nanstd(circadian_times_am),
                   np.nanstd(zeitgeber_times_pm),
                   np.nanstd(circadian_times_pm)],
             color=[utils.AM_COLOR, utils.AM_COLOR, utils.PM_COLOR, utils.PM_COLOR])

    axes.set_ylabel("Standard deviation (hours)", fontsize=14)

    x_ticks = [0, 1, 2, 3]
    x_tick_labels = ['Wall Clock', 'Biological Time',
                     'Wall Clock', 'Biological Time']

    for ax in [axes]:
        ax.set_ylim([0, y_bar_height + 0.1])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_tick_labels, fontsize=14)

    plt.tight_layout()
    plt.savefig("figures/compliance.png", dpi=utils.DPI)
    plt.show()


if __name__ == '__main__':
    with open(utils.snapshot_path, 'rb') as input_file:
        people = pickle.load(input_file)
    make_compliance_figure(people)
