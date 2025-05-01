import os
import pandas as pd
import AppKit
import time

import utils


def save_cleaned_data(subject, sorted_time, sorted_data):

    # Save base CSV
    df = pd.DataFrame({'Time': sorted_time, 'Data': sorted_data})
    csv_file_path = "data/cleaned/" + str(subject) + ".csv"
    df.to_csv(csv_file_path, index=False)

    # Save alternative format for Maria
    df['DaysSinceStart'] = (df['Time'] - df['Time'].min()) // (24 * 3600)
    df['HourOfDay'] = (df['Time'] % (24 * 3600)) // 3600
    df['Minutes'] = ((df['Time'] % (24 * 3600)) % 3600) // 60
    df = df[['DaysSinceStart', 'HourOfDay', 'Minutes', 'Data']]
    csv_file_path = "data/cleaned/" + str(subject) + "_v2.csv"
    df.to_csv(csv_file_path, index=False)


def raw_data_to_zeitgeber(data_frame, zeitgeber):
    extended_zeitgeber = []
    if utils.MODE == "light":
        extended_zeitgeber = zeitgeber + data_frame["LIGHT"].values.tolist()

    if utils.MODE == "activity":
        extended_zeitgeber = zeitgeber + \
            (utils.LIGHT_SCALAR * data_frame["PIM"].values).tolist()

    if utils.MODE == "light+activity":
        light_activity_hybrid = (
            utils.LIGHT_SCALAR * data_frame["PIM"].values) + data_frame["LIGHT"].values
        extended_zeitgeber = zeitgeber + light_activity_hybrid.tolist()

    return extended_zeitgeber


def load_and_clean_from_raw(subject):
    subject_directory = "data/Raw/" + subject + "/"
    all_epochs = []
    all_zeitgebers = []

    for file_name in os.listdir(subject_directory):
        file_with_path = os.path.join(subject_directory, file_name)

        if os.path.isfile(file_with_path) and file_with_path[-4:] == ".txt":
            print("Reading in " + str(file_with_path) + "...")

            with open(file_with_path, 'r') as file:
                first_line_check = file.readline()
            if "DATE" in first_line_check:
                data_frame = pd.read_csv(file_with_path, sep='\t')
            else:
                data_frame = pd.read_csv(file_with_path, sep='\t', skiprows=25)
            if 'DATE' not in data_frame.columns:
                data_frame = pd.read_csv(file_with_path, sep='\t', skiprows=24)
            data_frame["DateString"] = data_frame["DATE"] + \
                " " + data_frame["TIME"]

            data_frame['DateTime'] = pd.to_datetime(
                data_frame['DateString'], format="%d/%m/%Y %H:%M:%S", dayfirst=True, errors='coerce')
            data_frame['DateTime'].fillna(pd.to_datetime(
                data_frame['DateString'], format="%d/%m/%y %H:%M:%S", dayfirst=True, errors='coerce'), inplace=True)

            # Note: The below gives a "non-localized" error because there are hours in the dataset that
            # shouldn't exist in the America/Chicago timezone! (e.g. due to DST)
            # df['DateTime'] = df['DateTime'].dt.tz_localize('America/Chicago')

            first_day = data_frame['DateTime'].iloc[0].tz_localize(
                'UTC').tz_convert('America/Chicago')

            dst_offset = first_day.dst().total_seconds() if first_day.dst() else 0

            epoch_time = data_frame['DateTime'].view(
                'int64') // 1e9 - dst_offset
            all_epochs = all_epochs + epoch_time.tolist()
            all_zeitgebers = raw_data_to_zeitgeber(data_frame, all_zeitgebers)

    sorted_time = sorted(all_epochs)
    sorted_zeitgeber = [x for _, x in sorted(zip(all_epochs, all_zeitgebers))]
    save_cleaned_data(subject, sorted_time, sorted_zeitgeber)


if __name__ == '__main__':
    all_patients = ['P33', 'P34', 'P35', 'P36',
                    'P38', 'P39', 'P41', 'P42', 'P43', 'P44']
    start_time = time.time()

    for patient in all_patients:
        print(f"\nCleaning patient {patient}...")
        load_and_clean_from_raw(patient)

    execution_time = (time.time() - start_time)
    print('Execution time in seconds: ' + str(execution_time))
    AppKit.NSBeep()
