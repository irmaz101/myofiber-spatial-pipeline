"""
Annotate Cellpose training patches in Napari.

The script loads previously extracted image patches, allows manual
segmentation in Napari, and saves the resulting label masks. It can
also display patch locations on the original morphology image for
quality control.
"""

from pathlib import Path

import napari
import numpy as np
from skimage.io import imread, imsave


def load_patch(path):
    """Load an image if it exists."""

    path = Path(path)

    if path.exists():
        return imread(path)

    return None


def parse_patch_name(filepath):
    """Extract sample name and patch coordinates from the filename."""

    filename = Path(filepath).name

    sample_name = filename.split("-x")[0]

    xcoord = filename.split("-x")[1].split("-")[0].split("_")
    ycoord = filename.split("-y")[1].split(".")[0].split("_")

    xcoords = [int(x) for x in xcoord]
    ycoords = [int(y) for y in ycoord]

    square = np.array([
        [xcoords[0], ycoords[0]],
        [xcoords[0], ycoords[1]],
        [xcoords[1], ycoords[1]],
        [xcoords[1], ycoords[0]],
    ])

    return sample_name, square


def get_patch_location(patch_path):
    """Load the original image and return the annotated patch location."""

    sample_name, square = parse_patch_name(patch_path)

    image_path = Path("raw_images") / f"{sample_name}.tif"

    image = imread(image_path)

    return image, square


def show_patch_locations(image, squares):
    """Display patch locations on the original image."""

    viewer = napari.Viewer()

    viewer.add_image(image)

    viewer.add_shapes(
        squares,
        face_color="blue",
        edge_color="green",
        edge_width=3,
        name="training patches",
    )

    napari.run()


def annotate_patch(patch_path, label_path):
    """Open a patch in Napari and save the manually annotated mask."""

    patch = load_patch(patch_path)

    label = load_patch(label_path)

    if label is None:
        label = np.zeros(
            patch.shape[:2],
            dtype=np.uint16
        )

    viewer = napari.Viewer()

    viewer.add_image(patch)

    labels = viewer.add_labels(label)

    labels.brush_size = 4
    labels.mode = "paint"
    labels.preserve_labels = True

    napari.run()

    imsave(
        label_path,
        labels.data,
        check_contrast=False
    )


if __name__ == "__main__":
    patch_dir = Path("training_patches")
    label_dir = Path("training_labels")
    sample_pattern = "*"
    annotate = True

    squares = []

    for patch_path in sorted(patch_dir.glob(f"{sample_pattern}.tif")):

        label_path = (
            label_dir /
            f"{patch_path.stem}_label.png"
        )

        image, square = get_patch_location(patch_path)

        squares.append(square)

        if annotate:

            label_dir.mkdir(exist_ok=True)

            annotate_patch(
                patch_path,
                label_path
            )

    if not annotate and squares:

        show_patch_locations(
            image,
            squares
        )

    elif not squares:

        print("No patches found.")
