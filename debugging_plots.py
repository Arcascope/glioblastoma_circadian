import pickle
import matplotlib.pyplot as plt

import utils
from utils import *


def plot_std_over_time(participant_array, window=30):
    plt.close()
    for person in participant_array:
        color = [0.5, 0.5, 0.5]
        windowed = person.moving_std_with_window(window)
        plt.plot(windowed, 'o', alpha=0.4, color=color)
    plt.xlabel("Time (day)")
    plt.ylabel("SD in predicted DLMO")
    plt.savefig(output_folder + "dlmo_std_by_subject_" +
                str(window) + ".png", dpi=200)
    plt.close()


def plot_std_for_range_of_windows(participant_array):
    for window in [1, 5, 7, 14, 30, 60, 120, 180]:
        plot_std_over_time(participant_array=participant_array, window=window)


def generate_table_of_metrics(patient_array):
    print("\n\n Last four weeks DLMO standard deviation")
    for subject in patient_array:
        print(subject.subject_id + " " +
              str(subject.average_of_tail_with_length(4 * 7)))

    print("\n\nLast week DLMO standard deviation")
    for subject in patient_array:
        print(subject.subject_id + " " +
              str(subject.average_of_tail_with_length(1 * 7)))

    print("\n\nOverall DLMO standard deviation")
    for subject in patient_array:
        print(subject.subject_id + " " + str(subject.std_melatonin_onset))

    print("\n\nOverall pill timing standard deviation")
    for subject in patient_array:
        print(subject.subject_id + " " + str(subject.pill_std))

    print("\n\nOverall pill offset standard deviation")
    for subject in patient_array:
        print(subject.subject_id + " " + str(subject.pill_offset_std))

    with open(f"snapshots/{utils.MODE}.csv", 'w') as f:
        print("Subject, Last month DLMO STD, Last week DLMO STD, Overall DLMO STD, Pill STD, Pill Offset STD",
              file=f)
        for subject in patient_array:
            print(subject.subject_id + "," + str(subject.average_of_tail_with_length(4 * 7)) + "," + str(
                subject.average_of_tail_with_length(1 * 7)) + "," + str(subject.std_melatonin_onset) + "," + str(
                subject.pill_std) + "," + str(subject.pill_offset_std), file=f)


def plot_pill_melatonin_offset_group(participant_array):
    all_pill_offsets = []
    for subject in participant_array:
        all_pill_offsets += subject.all_pill_offsets.tolist()

    all_pill_offsets = np.array(all_pill_offsets)
    all_pill_offsets[all_pill_offsets <
                     5] = all_pill_offsets[all_pill_offsets < 5] + 24
    plt.close()

    _, _, _ = plt.hist(x=all_pill_offsets,
                       bins=np.arange(0, 36, 0.5),
                       color='#0504aa',
                       alpha=0.7,
                       rwidth=0.85)

    plt.savefig(output_folder + "overall_pill_offset_histogram.png", dpi=400)

    plt.close()


def plot_local_pill_timings_group(participant_array):
    all_pills_local_time = []
    for subject in participant_array:
        all_pills_local_time += subject.all_pills_local_time.tolist()

    all_pills_local_time = np.array(all_pills_local_time)
    all_pills_local_time[all_pills_local_time <
                         5] = all_pills_local_time[all_pills_local_time < 5] + 24
    plt.close()
    _, _, _ = plt.hist(x=all_pills_local_time,
                       bins=np.arange(0, 36, 0.5),
                       color='#b5040d',
                       alpha=0.7,
                       rwidth=0.85)
    plt.savefig(output_folder + "overall_pill_raw_histogram.png", dpi=400)


if __name__ == '__main__':
    with open(snapshot_path, 'rb') as input_file:
        people = pickle.load(input_file)

    generate_table_of_metrics(people)
    plot_pill_melatonin_offset_group(people)
    plot_local_pill_timings_group(people)
    plot_std_for_range_of_windows(people)
