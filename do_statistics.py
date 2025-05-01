import numpy as np
import pandas as pd
from scipy.stats import levene
import utils
from scipy.stats import kruskal
import pickle


def do_statistics(participant_array):
    zt_cut_point = 15
    ct_cut_point = 8

    zeitgeber_times_am = []
    zeitgeber_times_pm = []
    circadian_times_am = []
    circadian_times_pm = []

    zeitgeber_times_am_std = []  # Hold stds by person
    zeitgeber_times_pm_std = []
    circadian_times_am_std = []
    circadian_times_pm_std = []

    rows = []
    for subject in participant_array:
        subject_pills_local_time = subject.all_pills_local_time

        subject_pills_local_time[subject_pills_local_time > zt_cut_point] -= 24
        subject_pills_circadian_time = subject.all_pill_offsets
        subject_pills_circadian_time[subject_pills_circadian_time >
                                     ct_cut_point] -= 24

        subject_id = subject.subject_id
        rows.append([subject_id, "wall clock", subject.group] + [np.std(subject_pills_local_time)
                                                                 ] + [" ".join(map(str, subject_pills_local_time.tolist()))])
        rows.append([subject_id, "biological", subject.group] +
                    [np.std(subject_pills_circadian_time)] + [" ".join(map(str, subject_pills_circadian_time.tolist()))])

        if subject.group == utils.Group.AM:
            zeitgeber_times_am += subject_pills_local_time.tolist()
            zeitgeber_times_am_std += [np.std(subject_pills_local_time)]
            circadian_times_am += subject_pills_circadian_time.tolist()
            circadian_times_am_std += [np.std(
                subject_pills_circadian_time)]
        else:
            zeitgeber_times_pm += subject_pills_local_time.tolist()
            zeitgeber_times_pm_std += [np.std(subject_pills_local_time)]
            circadian_times_pm += subject_pills_circadian_time.tolist()
            circadian_times_pm_std += [np.std(subject_pills_circadian_time)]

    zeitgeber_times_am = np.array(zeitgeber_times_am)
    zeitgeber_times_pm = np.array(zeitgeber_times_pm)
    zeitgeber_times_am_std = np.array(zeitgeber_times_am_std)
    zeitgeber_times_pm_std = np.array(zeitgeber_times_pm_std)

    circadian_times_am = np.array(circadian_times_am)
    circadian_times_pm = np.array(circadian_times_pm)
    circadian_times_am_std = np.array(circadian_times_am_std)
    circadian_times_pm_std = np.array(circadian_times_pm_std)

    # Convert to DataFrame and round numbers to two decimal places
    df = pd.DataFrame(rows, columns=["Subject ID", "Time Type", "Group", "STD"] + [
                      f"Measurement {i+1}" for i in range(len(rows[0]) - 4)])
    df = df.round(2)
    # Print DataFrame to screen
    print(df.to_csv(index=False))

    # Write DataFrame to a CSV file
    output_file = "statistics_output.csv"
    df.to_csv(output_file, index=False)

    circadian_times_am = -circadian_times_am
    circadian_times_pm = -circadian_times_pm

    # Separate data into relevant groups
    am_wallclock = zeitgeber_times_am
    am_biological = circadian_times_am
    pm_wallclock = zeitgeber_times_pm
    pm_biological = circadian_times_pm

    # -------- 1. Levene's Test --------
    print("Levene's Test Results:")
    _, p_am = levene(am_wallclock, am_biological)
    _, p_pm = levene(pm_wallclock, pm_biological)
    print(f"AM Group: p-value = {p_am}")
    print(f"PM Group: p-value = {p_pm}")

    _, p_bio = levene(am_biological, pm_biological)
    _, p_wall = levene(am_wallclock, pm_wallclock)
    print(f"Biological Group: p-value = {p_bio}")
    print(f"Wall Clock Group: p-value = {p_wall}")

    # Compare across all groups
    _, p_all = levene(am_wallclock, am_biological, pm_wallclock, pm_biological)
    print(f"All Groups Comparison: p-value = {p_all}")

    # -------- 2. Bootstrapping Standard Deviations --------

    def bootstrap_std(data1, data2, n_bootstrap=1000):
        observed_diff = np.std(data1, ddof=1) - np.std(data2, ddof=1)
        combined = np.concatenate([data1, data2])
        bootstrap_diffs = []
        for _ in range(n_bootstrap):
            np.random.shuffle(combined)
            sample1 = combined[:len(data1)]
            sample2 = combined[len(data1):]
            bootstrap_diffs.append(
                np.std(sample1, ddof=1) - np.std(sample2, ddof=1))
        p_value = np.mean(np.abs(bootstrap_diffs) >= np.abs(observed_diff))
        return observed_diff, p_value

    print("\nBootstrapping Results:")
    obs_diff_am, p_am_boot = bootstrap_std(am_wallclock, am_biological)
    obs_diff_pm, p_pm_boot = bootstrap_std(pm_wallclock, pm_biological)
    print(
        f"AM Group: Observed Diff = {obs_diff_am:.3f}, p-value = {p_am_boot}")
    print(
        f"PM Group: Observed Diff = {obs_diff_pm:.3f}, p-value = {p_pm_boot}")

    print("\nKruskal-Wallis Test Results:")
    _, p_kruskal_am = kruskal(zeitgeber_times_am_std, circadian_times_am_std)
    _, p_kruskal_pm = kruskal(zeitgeber_times_pm_std, circadian_times_pm_std)
    print(f"AM Group: p-value = {p_kruskal_am}")
    print(f"PM Group: p-value = {p_kruskal_pm}")

    _, p_kruskal_bio = kruskal(circadian_times_am_std, circadian_times_pm_std)
    _, p_kruskal_wct = kruskal(zeitgeber_times_am_std, zeitgeber_times_pm_std)
    print(f"Bio Group: p-value = {p_kruskal_bio}")
    print(f"WCT Group: p-value = {p_kruskal_wct}")

    # Compare across all groups
    _, p_kruskal_all = kruskal(
        zeitgeber_times_am_std, circadian_times_am_std, zeitgeber_times_pm_std, circadian_times_pm_std)
    print(f"All Groups Comparison: p-value = {p_kruskal_all}")


if __name__ == "__main__":
    with open(utils.snapshot_path, 'rb') as input_file:
        people = pickle.load(input_file)
    do_statistics(people)
