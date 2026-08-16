"""
Extract random training patches from Xenium morphology images.

For each input image, the script extracts a specified number of random
1800 × 1800 pixel patches for manual annotation and Cellpose training.
"""

from pathlib import Path
from random import randint, seed

import numpy as np
from skimage.io import imread, imsave


def generate_random_patch(image, patch_size):
    """Generate a random square patch from an image."""

    height, width = image.shape[:2]

    if height < patch_size or width < patch_size:
        raise ValueError(
            f"Image shape {image.shape[:2]} is smaller than "
            f"the requested patch size {patch_size}."
        )

    x = randint(0, width - patch_size)
    y = randint(0, height - patch_size)

    if image.ndim == 3:
        patch = image[
            y:y + patch_size,
            x:x + patch_size,
            :
        ]
    else:
        patch = image[
            y:y + patch_size,
            x:x + patch_size
        ]

    idx_string = (
        f"-x{x}_{x + patch_size}"
        f"-y{y}_{y + patch_size}"
    )

    return patch, idx_string


def create_random_patch(image_path, patch_size):
    """Read an image and extract one random patch."""

    image = imread(image_path)
    return generate_random_patch(image, patch_size)


def save_patch(patch, sample_name, output_dir, idx_string):
    """Save a patch."""

    imsave(
        output_dir / f"{sample_name}{idx_string}.tif",
        patch,
        check_contrast=False
    )


if __name__ == "__main__":

    input_dir = Path("raw_images")
    output_dir = Path("training_patches")
    output_dir.mkdir(exist_ok=True)

    patches_per_image = 12
    patch_size = 1800
    random_seed = 0

    seed(random_seed)

    for image_path in sorted(input_dir.glob("*.tif")):

        sample_name = image_path.stem

        for _ in range(patches_per_image):

            patch, idx_string = create_random_patch(
                image_path,
                patch_size
            )

            save_patch(
                patch,
                sample_name,
                output_dir,
                idx_string
            )