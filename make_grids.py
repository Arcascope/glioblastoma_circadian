import os
from PIL import Image, ImageDraw, ImageFont
from utils import subject_id_to_number


def make_grids_for_ext(ext, draw_label=False):
    all_patients = ['P33', 'P34', 'P35', 'P36',
                    'P38', 'P39', 'P41', 'P42', 'P43', 'P44']

    images = []
    for patient in all_patients:
        image_path = os.path.join('debug_outputs/', f'{patient}{ext}.png')
        if os.path.exists(image_path):
            img = Image.open(image_path)
            if draw_label:
                padding = 150
                img_with_padding = Image.new(
                    'RGB', (img.size[0], img.size[1] + padding), 'white')
                img_with_padding.paste(img, (0, padding))
                images.append(img_with_padding)
            else:
                images.append(img)
        else:
            print(f"Image not found: {image_path}")
    if images:
        grid_width = 3
        grid_height = 4

        img_width = max(img.size[0] for img in images)
        img_height = max(img.size[1] for img in images)

        combined_image = Image.new(
            'RGB', (img_width * grid_width, img_height * grid_height), 'white')

        for i, img in enumerate(images):

            x_offset = (i % grid_width) * img_width
            y_offset = (i // grid_width) * img_height
            combined_image.paste(images[i], (x_offset, y_offset))

            if draw_label:
                label = f'P{subject_id_to_number[all_patients[i]]}'
                draw = ImageDraw.Draw(combined_image)
                font = ImageFont.truetype("/Library/Fonts/Arial.ttf", 200)
                draw.text((x_offset + 10, y_offset),
                          label, fill='black', font=font)

        combined_image.save(f'debug_outputs/combined_grid{ext}.png')


def make_grids():
    make_grids_for_ext("", draw_label=True)
    make_grids_for_ext("_pill_histogram")
    make_grids_for_ext("_melatonin_onset_histogram")
    make_grids_for_ext("_pill_offset_histogram")
    make_grids_for_ext("_scatter", draw_label=True)


if __name__ == '__main__':
    make_grids()
