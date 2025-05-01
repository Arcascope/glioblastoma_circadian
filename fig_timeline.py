import pandas as pd
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import utils
import numpy as np


def make_timeline_figure():
    forward_hatch = '///////'
    backward_hatch = '\\\\\\\\\\\\\\'
    # Load the data
    data = utils.read_outcomes_spreadsheet()
    event_colors = {
        'Date of Pseudoprogression': utils.GOLD_REGULAR,
        'Date of Progression (1)': utils.ORANGE_REGULAR,
        'Date of Progression (2)': utils.RED_REGULAR
    }

    # Set up the plot
    # Adjust the size as necessary
    fig, ax = plt.subplots(figsize=(10, len(data) * 0.45))
    # Customize the plot
    ax.set_yticks(np.array(range(len(data))) + 0.4)

    ids = [utils.subject_id_to_number["P" + str(id_value)]
           for id_value in data['ID'].values]
    ax.set_yticklabels(ids, fontsize=13)
    ax.set_xlabel('Days since Surgery', fontsize=13)
    ax.set_xlim([0, 2200])
    ax.set_ylim([0, len(data)])

    y = 0  # Initial y position
    row_height = 0.8  # Adjust as needed based on your plot

    # Draw the rectangles and markers
    for _, row in data.iterrows():
        # Lifespan rectangle
        lifespan_width = (row['Date of Death'] - row['Surgery date']).days
        print(lifespan_width)

        if row['Arm'] == "AM":
            hatch = forward_hatch
            drug_base_color = utils.AM_COLOR
        else:
            hatch = backward_hatch
            drug_base_color = utils.PM_COLOR

        ax.add_patch(
            patches.Rectangle((0, y), lifespan_width, 0.8, color=utils.BACKGROUND_COLOR,
                              edgecolor=[0, 0, 0, 0]))

        if row['Alive']:
            # Adding a small offset to the right of the rectangle
            text_x_position = 0 + lifespan_width + 20
            text_y_position = y + 0.35  # Vertically centered on the rectangle
            ax.text(text_x_position, text_y_position, "A",
                    va='center', ha='left', fontsize=20)

        # Drug duration rectangle
        drug_start = (row['Start Date of Drug'] - row['Surgery date']).days
        drug_end = (row['End Date of Drug'] - row['Surgery date']).days
        ax.add_patch(
            patches.Rectangle((drug_start, y), drug_end - drug_start, 0.8, fill=None, hatch=hatch,
                              color=drug_base_color,
                              edgecolor=tuple(
                                  element * 0 for element in drug_base_color)
                              ))

        # Markers for key events with different colors and ensuring circular shape
        for event, color in event_colors.items():
            if pd.notnull(row[event]):
                event_day = (row[event] - row['Surgery date']).days
                ellipse = patches.Ellipse(
                    (event_day, y + 0.4), width=60, height=row_height * 0.74, color=color)
                ax.add_patch(ellipse)
                ellipse = patches.Ellipse((event_day, y + 0.4), width=60, height=row_height * 0.74,
                                          color=[0, 0, 0, 0.25],
                                          fill=False)
                ax.add_patch(ellipse)

        y += 1  # Move to the next row

    legend_elements = [Line2D([0], [0], marker='o', color='w', label=event,
                              markersize=10, markerfacecolor=color) for event, color in event_colors.items()]

    # Add legend elements for hatched rectangles
    legend_elements += [patches.Patch(facecolor=utils.BACKGROUND_COLOR, edgecolor=utils.AM_COLOR,
                                      label="AM Dosing Period", hatch=forward_hatch)]
    legend_elements += [patches.Patch(facecolor=utils.BACKGROUND_COLOR, edgecolor=utils.PM_COLOR,
                                      label="PM Dosing Period", hatch=backward_hatch)]

    ax.legend(handles=legend_elements)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.yaxis.set_tick_params(length=0)
    plt.savefig("figures/timeline.png", dpi=utils.DPI)
    plt.show()


if __name__ == '__main__':
    make_timeline_figure()
