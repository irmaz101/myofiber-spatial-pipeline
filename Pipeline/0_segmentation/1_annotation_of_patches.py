"""Annotate training patches in Napari, or inspect their image locations."""

from collections import defaultdict
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
    """Extract sample name and patch coordinates from the filename.

    The returned square uses the Napari (row, column) convention,
    that is (y, x).
    """

    filename = Path(filepath).name

    sample_name = filename.split("-x")[0]

    xcoord = filename.split("-x")[1].split("-")[0].split("_")
    ycoord = filename.split("-y")[1].split(".")[0].split("_")

    xcoords = [int(x) for x in xcoord]
    ycoords = [int(y) for y in ycoord]

    square = np.array(
        [
            [ycoords[0], xcoords[0]],
            [ycoords[1], xcoords[0]],
            [ycoords[1], xcoords[1]],
            [ycoords[0], xcoords[1]],
        ]
    )

    return sample_name, square


def show_patch_locations(image, squares, sample_name):
    """Display patch locations on the original image."""

    viewer = napari.Viewer()

    viewer.add_image(image, name=sample_name)

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
        label = np.zeros(patch.shape[:2], dtype=np.uint16)

    viewer = napari.Viewer()

    viewer.add_image(patch)

    labels = viewer.add_labels(label)

    labels.brush_size = 4
    labels.mode = "paint"
    labels.preserve_labels = True

    napari.run()

    imsave(label_path, labels.data, check_contrast=False)


if __name__ == "__main__":

    image_dir = Path("raw_images")
    patch_dir = Path("training_patches")
    label_dir = Path("training_labels")
    sample_pattern = "*"

    # True: annotate each patch in Napari.
    # False: show the patch locations on the original images.
    annotate = True

    patch_paths = sorted(patch_dir.glob(f"{sample_pattern}.tif"))

    if not patch_paths:
        print("No patches found.")

    elif annotate:

        label_dir.mkdir(exist_ok=True)

        for patch_path in patch_paths:

            annotate_patch(patch_path, label_dir / f"{patch_path.stem}_label.png")

    else:

        squares_by_sample = defaultdict(list)

        for patch_path in patch_paths:
            sample_name, square = parse_patch_name(patch_path)
            squares_by_sample[sample_name].append(square)

        for sample_name, squares in squares_by_sample.items():

            show_patch_locations(
                imread(image_dir / f"{sample_name}.tif"), squares, sample_name
            )
