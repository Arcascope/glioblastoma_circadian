import pickle
from scipy.optimize import minimize

import matplotlib.pyplot as plt
import utils
from PIL import Image, ImageDraw, ImageFont
import numpy as np
from matplotlib import font_manager


def make_scatter_plot_figure(participant_array, save_name="scatter"):
    zt_cut_point = 15
    ct_cut_point = 8
    has_legend = {"AM": False, "PM": False}
    all_circadian_times = np.array([])
    all_zeitgeber_times = np.array([])

    # Scatter plot: ZT vs CT for all pill timings.
    for subject in participant_array:
        zeitgeber_time = subject.all_pills_local_time  # Local time of pill intake
        # Hours before DLMO pill timing
        circadian_time = subject.all_pill_offsets

        zeitgeber_time[zeitgeber_time > zt_cut_point] -= 24
        circadian_time[circadian_time > ct_cut_point] -= 24

        # Change to hours after DLMO, so 1 hour before DLMO becomes 23 hours after DLMO
        circadian_time = 24 - circadian_time

        all_zeitgeber_times = np.append(all_zeitgeber_times, zeitgeber_time)
        all_circadian_times = np.append(all_circadian_times, circadian_time)
        plt.scatter(zeitgeber_time, circadian_time, marker='o', color=subject.group.color(),
                    label=f"{subject.group.value} Dose" if not has_legend[subject.group.value] else None)

        has_legend[subject.group.value] = True

    def objective(c):
        # Calculate residuals with fixed slope of 1
        residuals = all_circadian_times - (all_zeitgeber_times + c)
        return np.sum(residuals ** 2)

    # Minimize the objective function to find the best intercept c
    result = minimize(objective, x0=0)  # Start with an initial guess for c
    best_intercept = result.x[0]
    x_points = np.array([np.min(all_zeitgeber_times) - 2,
                        np.max(all_zeitgeber_times) + 2])

    # Calculate RMS error for the fit
    residuals = all_circadian_times - (all_zeitgeber_times + best_intercept)
    rms_error = np.sqrt(np.mean(residuals ** 2))
    max_error = np.max(np.abs(residuals))
    print(f"RMS Error: {rms_error}")
    print(f"Max Error: {max_error}")

    # Calculate the number of points more than 2 hours away from the line of best fit
    threshold = 1
    num_points_far_away = np.sum(np.abs(residuals) > threshold)
    print(
        f"Number of points more than {threshold} hours away from the line of best fit: {num_points_far_away}")
    total_points = len(all_circadian_times)
    percent_far_away = (num_points_far_away / total_points) * 100
    print(
        f"Percentage of points more than {threshold} hours away from the line of best fit: {percent_far_away:.2f}%")

    # Calculate R^2 value
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((all_circadian_times - np.mean(all_circadian_times)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    print(f"R^2: {r_squared}")

    plt.plot(x_points, best_intercept + x_points, 'k--')

    label_font_size = 14
    tick_font_size = 11
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    plt.xlabel('Wall clock time', fontsize=label_font_size)
    plt.ylabel('Biological time (hours after DLMO mod 24)',
               fontsize=label_font_size)

    plt.xticks(np.arange(-8, 16, 4))
    current_xticks = plt.xticks()[0]
    hour_labels = ["{:02d}:00".format(
        round(np.mod(hour + 24, 24)))for hour in current_xticks]

    plt.xticks(current_xticks, hour_labels, fontsize=tick_font_size)

    current_yticks = plt.yticks()[0]
    y_labels = [np.mod(hour, 24)
                for hour in current_yticks]
    plt.yticks(current_yticks, y_labels, fontsize=tick_font_size)

    plt.text(0.65, 0.22, f'RMS Error: {rms_error:.1f} hrs', transform=plt.gca().transAxes,
             fontsize=16, verticalalignment='top')
    plt.text(0.65, 0.15, f'Max Error: {max_error:.1f} hrs', transform=plt.gca().transAxes,
             fontsize=16, verticalalignment='top')
    plt.text(0.65, 0.08, f'$R^2$: {r_squared:.2f}', transform=plt.gca().transAxes,
             fontsize=16, verticalalignment='top')

    plt.yticks(fontsize=tick_font_size)
    plt.legend()
    if save_name == "scatter":
        plt.savefig(f"figures/{save_name}.png", dpi=utils.DPI)
    else:
        plt.savefig(f"debug_outputs/{save_name}_scatter.png", dpi=utils.DPI)
    plt.show()


def pad_image(original_image, padding_width=None, padding_height=None):
    # Calculate the new image size
    if padding_width is not None:
        new_width = original_image.width + 2 * padding_width
    else:
        new_width = original_image.width

    if padding_height is not None:
        new_height = original_image.height + 2 * padding_height
    else:
        new_height = original_image.height

    # Create a new image with white padding
    padded_image = Image.new('RGB', (new_width, new_height), 'white')

    # Paste the original image onto the padded image with padding offsets
    x_offset = padding_width if padding_width is not None else 0
    y_offset = padding_height if padding_height is not None else 0
    padded_image.paste(original_image, (x_offset, y_offset))
    return padded_image


def combine_actogram_and_scatterplot():
    image1 = Image.open('figures/scatter.png')
    image2 = Image.open('debug_outputs/P39.png')

    image1 = pad_image(image1, padding_height=100)

    # Get the dimensions of both images
    width1, height1 = image1.size
    width2, height2 = image2.size

    # Determine the taller image
    if height1 > height2:
        taller_image = image1
        shorter_image = image2
    else:
        taller_image = image2
        shorter_image = image1

    # Calculate the scaling factor to match the height of the taller image
    scaling_factor = taller_image.height / shorter_image.height
    new_width = int(shorter_image.width * scaling_factor)

    # Resize the shorter image to match the height of the taller image while preserving aspect ratio
    resized_shorter_image = shorter_image.resize(
        (new_width, taller_image.height))

    # Create a new image with white padding on the left
    combined_width = taller_image.width + resized_shorter_image.width
    combined_height = taller_image.height
    combined_image = Image.new(
        'RGB', (combined_width, combined_height), 'white')

    # Paste the resized images onto the combined image
    combined_image.paste(resized_shorter_image, (0, 0))
    combined_image.paste(taller_image, (resized_shorter_image.width, 0))

    draw = ImageDraw.Draw(combined_image)
    font_size = utils.DPI / 2.25
    buffer = font_size

    # Find the path for a specific font:
    file = font_manager.findfont('Helvetica')

    # You can change the font family if needed
    font = ImageFont.truetype(file, font_size)
    draw.text((buffer / 2, buffer), "A", fill='black', font=font)
    draw.text((resized_shorter_image.width + buffer /
              2, buffer), "B", fill='black', font=font)

    combined_image.save('figures/combined_scatter_plot_and_actogram.png')
    # combined_image.show()


if __name__ == '__main__':
    with open(utils.snapshot_path, 'rb') as input_file:
        people = pickle.load(input_file)
    make_scatter_plot_figure(people)
    combine_actogram_and_scatterplot()
