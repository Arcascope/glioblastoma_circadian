import pickle
from datetime import datetime
import matplotlib.pyplot as plt
import pytz
import time
import utils
from person_metrics import PersonMetrics
from simulation_result import SimulationResult
from utils import *
from lco import SinglePopModel
import matplotlib
from scipy.ndimage import gaussian_filter1d

pd.set_option('display.float_format', lambda x: '%.3f' % x)

DELTA_T = 0.1
DAYS_POST_MISSING_TO_IGNORE = 4
MINUTES_PER_HOUR = 60
SECONDS_PER_HOUR = 3600
HOURS_PER_DAY = 24
INVALID_DAY_THRESHOLD = 50
x_label_font_size = 18
title_font_size = 46
x_tick_font_size = 16
y_tick_font_size = 16
use_same_axes = True

chicago_timezone = pytz.timezone('America/Chicago')


# Return np.nan if all values are nan; else, return the nanmean.
def extended_nanmean(array_slice):
    if np.all(np.isnan(array_slice)):
        return np.nan
    else:
        return np.nanmean(array_slice)


def run_model(subject_id, timestamps, zeitgeber, run_mode=""):
    # Use UTC to avoid issues with DST (e.g. "non-existent dates") [America/Chicago does not work]
    tz = pytz.timezone('UTC')

    dt = datetime.fromtimestamp(timestamps[0], tz)
    hour = int(dt.strftime('%H')) + int(dt.strftime('%M')) / MINUTES_PER_HOUR

    initial_condition = np.array([0.6, phase_ic_guess(hour), 0.0])

    # Make zero the first timestamp and convert to hours
    timestamps = timestamps - timestamps[0]
    timestamps = timestamps / SECONDS_PER_HOUR

    model = SinglePopModel()

    # Zero out NaN light or activity levels
    zeitgeber[np.isnan(zeitgeber)] = 0

    # TODO: Replace with circadian package
    sol = model.integrate_model(timestamps,
                                zeitgeber,
                                initial_condition)

    plt.xlabel('Time hours since start')
    plt.ylabel('Model output')
    plt.plot(timestamps, np.multiply(sol[0, :], np.cos(sol[1, :])))
    plt.savefig(output_folder + str(subject_id) +
                "_" + run_mode + "_model_output.png")
    plt.close()

    melatonin_onset_times = model.integrate_observer(timestamps,
                                                     zeitgeber,
                                                     initial_condition,
                                                     SinglePopModel.DLMOObs)

    plt.plot(np.mod(hour + melatonin_onset_times, HOURS_PER_DAY), 'go')

    plt.savefig(output_folder + str(subject_id) + "_" +
                run_mode + "_debug_melatonin_plot.png")
    plt.close()
    return sol, melatonin_onset_times


def save_circadian_amplitude(time, sol, zeitgeber, dt, save_name, invalid_days):
    circadian_amplitude = sol[0, :]

    # Calculate the daily average of circadian amplitude
    days_with_data = np.arange(
        0, len(circadian_amplitude) * dt / SECONDS_PER_HOUR / HOURS_PER_DAY)

    # Pop last element, to match SRI calculations
    days_with_data = days_with_data[:-1]
    daily_circadian_amplitude = []

    for day in days_with_data:
        start_index = int(day * SECONDS_PER_HOUR * HOURS_PER_DAY / dt)
        end_index = int((day + 1) * SECONDS_PER_HOUR * HOURS_PER_DAY / dt)
        if end_index >= len(circadian_amplitude):
            end_index = len(circadian_amplitude) - 1
        if start_index != end_index:
            daily_avg = np.mean(circadian_amplitude[start_index:end_index])

            if np.isnan(np.sum(zeitgeber[start_index:end_index])):
                daily_circadian_amplitude.append(np.nan)
            else:
                daily_circadian_amplitude.append(daily_avg)
        else:
            daily_circadian_amplitude.append(np.nan)

    # Confirm remove all invalid days
    for day in invalid_days:
        if day in days_with_data:
            index = np.where(days_with_data == day)[0][0]
            daily_circadian_amplitude[index] = np.nan

    # Save circadian amplitude to a CSV file
    circadian_amplitude_df = pd.DataFrame({
        "Days with Data": days_with_data,
        "Circadian Amplitude": daily_circadian_amplitude
    })

    circadian_amplitude_df.to_csv(
        output_folder + save_name + "_circadian_amplitude.csv", index=False)

    # Plot daily circadian amplitude over time
    plt.plot(days_with_data, daily_circadian_amplitude,
             'bo', label='Daily Circadian Amplitude')
    plt.xlabel('Days with Data')
    plt.ylabel('Daily Circadian Amplitude')
    plt.title('Daily Circadian Amplitude Over Time')
    plt.legend()
    plt.savefig(output_folder + save_name +
                "_daily_circadian_amplitude_plot.png")
    plt.close()
    return daily_circadian_amplitude


def save_predictive_csvs(epochs, zeitgeber, dt, save_name, invalid_days, sol):
    _ = save_circadian_amplitude(epochs,
                                 sol,
                                 zeitgeber,
                                 dt,
                                 save_name,
                                 [])


def run_patient_and_plot(epochs, zeitgeber, save_name, local_pill_times=None, dst_corrected_pill_times=None):
    if local_pill_times is None or dst_corrected_pill_times is None:
        local_pill_times = []
        dst_corrected_pill_times = []

    min_time = np.min(epochs)
    max_time = np.max(epochs)
    epochs = np.array(epochs)
    zeitgeber = np.array(zeitgeber)
    dt = epochs[1] - epochs[0]

    # .copy() makes sure we don't run into any pass-by-reference issues.
    sol, melatonin_onset_times = run_model(
        save_name, epochs, zeitgeber.copy(), run_mode=utils.MODE)
    num_days = int(np.floor((max_time - min_time) /
                   (SECONDS_PER_HOUR * HOURS_PER_DAY)))

    # For all (DST-corrected) pill times, identify the distance from the nearest melatonin onset
    # If there is valid pill timing then, save both offset and local time.
    pill_offsets = []
    valid_local_pill_times = []
    for corrected_pill_time, local_pill_time in zip(dst_corrected_pill_times, local_pill_times):
        melatonin_onset_in_seconds = min_time + melatonin_onset_times * SECONDS_PER_HOUR
        index = np.argmin(
            np.abs(melatonin_onset_in_seconds - corrected_pill_time))

        # Offsets here is in units of "hours before DLMO"
        # One hour after melatonin onset would give us -1, etc.
        closest_difference = (
            melatonin_onset_in_seconds[index] - corrected_pill_time) / SECONDS_PER_HOUR

        # Only save pills with valid actigraphy data within a day
        if np.abs(closest_difference) <= 24:
            pill_offsets.append(np.mod(closest_difference + 24, 24))
            valid_local_pill_times.append(local_pill_time)

    print(f"\nHow many pills did we not have valid data for {save_name}?")
    print("Total pills: " + str(len(dst_corrected_pill_times)))
    print("Valid pills: " + str(len(valid_local_pill_times)))

    actogram = []
    invalid_days = []
    pills_on_days = []

    # Create an actogram
    for i in range(num_days):
        activity_in_day = []
        start_of_day = min_time + HOURS_PER_DAY * SECONDS_PER_HOUR * i
        end_of_day = min_time + HOURS_PER_DAY * SECONDS_PER_HOUR * (i + 1)

        # Get pill intakes within day
        pills_for_day = dst_corrected_pill_times[
            (dst_corrected_pill_times >= start_of_day) & (dst_corrected_pill_times < end_of_day)]

        if len(pills_for_day) > 0:
            pills = [np.mod(pill, HOURS_PER_DAY * SECONDS_PER_HOUR) /
                     SECONDS_PER_HOUR for pill in pills_for_day]
            pills_on_days.append(pills)
        else:
            pills_on_days.append([np.nan])

        # Iterate over chunks of day
        for time_chunk in range(int(HOURS_PER_DAY / DELTA_T)):
            time_of_interest = start_of_day + time_chunk * DELTA_T * SECONDS_PER_HOUR
            steps_in_range = zeitgeber[
                (epochs >= time_of_interest) & (epochs < time_of_interest + DELTA_T * SECONDS_PER_HOUR)]
            if len(steps_in_range) > 0:
                ac_value = extended_nanmean(steps_in_range)
            else:
                ac_value = np.nan
            activity_in_day.append(ac_value)

        # If a day has no data (or is below a threshold), mark it as invalid
        if np.sum(activity_in_day) < INVALID_DAY_THRESHOLD or np.isnan(np.sum(activity_in_day)):
            invalid_days += [i]
            actogram.append(np.ones_like(activity_in_day) * np.nan)
        else:
            actogram.append(activity_in_day)

    actogram_shifted = actogram[1:]
    actogram_shifted.append(actogram[0])
    actogram = np.hstack((np.array(actogram), np.array(actogram_shifted)))

    invalid_days += [0]

    # Burn in days we want to ignore after invalid days
    for i in range(DAYS_POST_MISSING_TO_IGNORE):
        invalid_days += [day + 1 for day in invalid_days]

    # Remove duplicates from previous step
    invalid_days = np.unique(invalid_days)

    # Ignore any invalid days that extend past the end of the actogram
    invalid_days = invalid_days[invalid_days < np.shape(actogram)[0]]
    melatonin_onset_times[invalid_days] = np.nan

    # Save data for analysis of progression
    save_predictive_csvs(epochs, zeitgeber, dt, save_name, invalid_days, sol)

    pill_offsets = np.array(pill_offsets)

    # Rescale and plot
    actogram = np.log(actogram + 1)
    cmap = matplotlib.colormaps["viridis"]
    cmap.set_bad(color='k')  # Show NaNs in different color
    _ = plt.figure(figsize=(10, 1.4 + 0.05 * len(actogram) / 2))

    plt.imshow(actogram, interpolation='none',
               aspect='auto', vmin=0, vmax=10, cmap=cmap)

    # Plot pill timing (corrected for DST)
    y_range = list(range(len(melatonin_onset_times)))
    x_range = np.mod(melatonin_onset_times - 12, HOURS_PER_DAY) + 12
    y_range_pills = list(range(len(pills_on_days)))
    plt.plot(x_range / DELTA_T, y_range, '--', color=utils.RED_REGULAR)

    for pill_list, y_location in zip(pills_on_days, y_range_pills):
        pill_list = np.array(pill_list)
        for pill_value in pill_list:
            plt.plot(pill_value / DELTA_T, y_location, 'o',
                     markerfacecolor=utils.ORANGE_LIGHT,
                     markeredgecolor=utils.ORANGE_REGULAR)

    plt.xticks([0, 12 / DELTA_T, 24 / DELTA_T, 36 / DELTA_T, 48 / DELTA_T],
               ["00:00", "12:00", "00:00", "12:00", "00:00"], fontsize=x_tick_font_size)
    plt.yticks(fontsize=y_tick_font_size)
    plt.xlabel("Standard time", fontsize=x_label_font_size)
    plt.ylabel("Days since study start", fontsize=x_label_font_size)

    plt.tight_layout(pad=2.0)

    plt.savefig(output_folder + save_name + ".png", dpi=utils.DPI)
    plt.close()

    # Convert raw pill times (no DST correction) to local time
    valid_local_pill_times = np.array(valid_local_pill_times)
    valid_local_pill_times = np.mod(
        valid_local_pill_times, HOURS_PER_DAY * SECONDS_PER_HOUR) / SECONDS_PER_HOUR

    plot_local_pill_time_histogram(valid_local_pill_times, save_name)
    plot_offset_histograms(melatonin_onset_times, pill_offsets, save_name)

    simulation_result = SimulationResult(sol=sol,
                                         dlmos=melatonin_onset_times,
                                         zeitgeber=zeitgeber,
                                         pills_local_time=valid_local_pill_times,
                                         pills_corrected_for_dst=dst_corrected_pill_times,
                                         pill_melatonin_offsets=pill_offsets,
                                         actogram=actogram)
    return simulation_result


def plot_local_pill_time_histogram(local_pill_times, save_name):
    _, _, _ = plt.hist(x=np.mod(np.mod(local_pill_times - 12, HOURS_PER_DAY) + 12, 24),
                       bins=np.arange(0, 24, 0.5),
                       color='#dd040a',
                       alpha=0.9,
                       rwidth=0.85)
    plt.xlabel("Hour of day (wall clock time)", fontsize=x_label_font_size)
    plt.title(
        "P"+str(utils.subject_id_to_number[save_name]), fontsize=title_font_size)
    plt.xticks(np.arange(0, 25, 2), fontsize=x_tick_font_size)
    plt.yticks(fontsize=y_tick_font_size)
    if use_same_axes:
        plt.ylim([0, 40])
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout(pad=2.0)
    plt.savefig(output_folder + save_name + "_pill_histogram.png", dpi=300)
    plt.close()


def plot_offset_histograms(melatonin_onset_times, pill_offsets, save_name):
    _, _, _ = plt.hist(x=np.mod(np.mod(melatonin_onset_times - 12, HOURS_PER_DAY) + 12, 24),
                       bins=np.arange(0, 24, 0.5),
                       color='#0504aa',
                       alpha=0.7,
                       rwidth=0.85)
    plt.xlabel("Hour of day (standard time)", fontsize=x_label_font_size)
    plt.title(
        "P"+str(utils.subject_id_to_number[save_name]), fontsize=title_font_size)
    plt.xticks(np.arange(0, 25, 2), fontsize=x_tick_font_size)
    plt.yticks(fontsize=y_tick_font_size)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    if use_same_axes:
        plt.ylim([0, 90])

    plt.tight_layout(pad=2.0)
    plt.savefig(output_folder + save_name +
                "_melatonin_onset_histogram.png", dpi=300)
    plt.close()

    _, _, _ = plt.hist(x=24 - np.mod(np.mod(pill_offsets - 12, HOURS_PER_DAY) + 12, 24),
                       bins=np.arange(0, 24, 0.5),
                       color='#e504a6',
                       alpha=0.9,
                       rwidth=0.85)
    plt.xlabel("Hours after DLMO", fontsize=x_label_font_size)
    plt.title(
        "P"+str(utils.subject_id_to_number[save_name]), fontsize=title_font_size)
    plt.xticks(np.arange(0, 25, 2), fontsize=x_tick_font_size)
    plt.yticks(fontsize=y_tick_font_size)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    if use_same_axes:
        plt.ylim([0, 20])

    plt.tight_layout(pad=2.0)
    plt.savefig(output_folder + save_name +
                "_pill_offset_histogram.png", dpi=300)


def prepend_midnight_and_run(subject, sorted_time, sorted_zeitgeber):
    start_date = datetime.utcfromtimestamp(np.min(sorted_time))
    midnight_on_start_day = (start_date - start_date.replace(hour=0,
                                                             minute=0,
                                                             second=0,
                                                             microsecond=0)).total_seconds()

    # Ensure all actograms start at midnight (Note: Assuming UTC for all times here)
    sorted_time = np.insert(sorted_time, 0, np.min(
        sorted_time) - midnight_on_start_day)
    sorted_zeitgeber = np.insert(sorted_zeitgeber, 0, 0)

    timestep_in_seconds = DELTA_T * SECONDS_PER_HOUR
    time_interp = np.arange(np.min(sorted_time), np.max(
        sorted_time), timestep_in_seconds)

    sorted_time = np.array(sorted_time)
    sorted_zeitgeber = np.array(sorted_zeitgeber)

    zeitgeber_interpolated = []
    progress_count = 0
    for time_slice in time_interp:
        progress_count = progress_count + 1
        if progress_count % 10000 == 0:
            print(str(100 * progress_count / len(time_interp)) + "% complete...")

        indices = np.where((sorted_time >= time_slice) & (
            sorted_time < time_slice + timestep_in_seconds))[0]
        zeitgeber_in_range = sorted_zeitgeber[indices]

        if len(zeitgeber_in_range) == 0:
            zeitgeber_interpolated.append(np.nan)
        else:
            zeitgeber_interpolated.append(extended_nanmean(zeitgeber_in_range))

    # See README for a note about DST in the below.
    pill_diaries = pd.read_excel(
        'data/PillDiaries.xlsx', sheet_name=str(subject[1:]))

    # Convert 'Date' and 'What time was dose taken' columns to datetime objects
    date_column = pd.to_datetime(pill_diaries['Date'])
    time_column = pd.to_datetime(
        pill_diaries['What time was dose taken'], format='%H:%M:%S').dt.time

    # Combine date and time components into a single datetime column
    combined_dates = date_column + pd.to_timedelta(time_column.astype(str))

    localized_dates = combined_dates.apply(
        lambda x: chicago_timezone.localize(x))
    pill_dst_correction = [dt.astimezone(
        chicago_timezone).dst().total_seconds() for dt in localized_dates]

    pill_times_local = (combined_dates.view('int64') // 1e9).values
    dst_corrected_pill_times = pill_times_local - pill_dst_correction
    simulation_result = run_patient_and_plot(time_interp,
                                             zeitgeber_interpolated,
                                             save_name=str(subject),
                                             local_pill_times=pill_times_local,
                                             dst_corrected_pill_times=dst_corrected_pill_times)
    return simulation_result


def save_output(obj, filename):
    with open(filename, 'wb') as output_file:  # Overwrites any existing file.
        pickle.dump(obj, output_file, pickle.HIGHEST_PROTOCOL)


def load_and_run_cleaned(subject):
    file_path = f"data/cleaned/{subject}.csv"
    df = pd.read_csv(file_path)

    sorted_time = df['Time'].to_numpy()
    sorted_zeitgeber = df['Data'].to_numpy()

    simulation_result = prepend_midnight_and_run(
        subject, sorted_time, sorted_zeitgeber)
    return simulation_result


if __name__ == '__main__':

    start_time = time.time()
    all_patients = ['P33', 'P34', 'P35', 'P36',
                    'P38', 'P39', 'P41', 'P42', 'P43', 'P44']
    # all_patients = ['P33']
    data = read_outcomes_spreadsheet()

    people = []
    for patient in all_patients:
        print(f"\nRunning patient {patient}...")

        # Run simulation
        simulation_for_person = load_and_run_cleaned(patient)

        # Find what group they were assigned to and store
        arm = data[data['ID'] == int(patient[1:])]['Arm'].iloc[0]
        group = Group.AM if arm == "AM" else Group.PM

        # Add to the list to save
        people.append(PersonMetrics(patient, group, simulation_for_person))

    save_output(people, snapshot_path)

    execution_time = (time.time() - start_time)
    print('Execution time in seconds: ' + str(execution_time))
