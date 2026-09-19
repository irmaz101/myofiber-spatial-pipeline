"""Segment full images with Cellpose-SAM and save masks at both image scales."""

from cellpose import models
from tifffile import imread, imwrite
from skimage.transform import resize
import numpy as np
import torch
from pathlib import Path
import gc

# Settings

input_dir = Path("arhiv_crop")
output_dir = Path("cellpose_output")
output_dir.mkdir(parents=True, exist_ok=True)

scale_factor = 0.5
diameter_orig = None

# Set to None to use the default pretrained Cellpose model.
# To use a custom model, provide the path to the trained model weights.
model_path = None
# model_path = Path("models/myofibre_model")

gpu_available = torch.cuda.is_available()


# Load model

if model_path is None:
    print("Using the default pretrained Cellpose model.")
    model = models.CellposeModel(gpu=gpu_available)
else:
    if not Path(model_path).exists():
        raise FileNotFoundError(f"Custom model not found: {model_path}")

    print(f"Using custom Cellpose model: {model_path}")
    model = models.CellposeModel(gpu=gpu_available, pretrained_model=str(model_path))


diameter_rescaled = None if diameter_orig is None else diameter_orig * scale_factor


# Segment images

for tif_path in sorted(input_dir.glob("*.tif")):
    print(f"\nProcessing {tif_path.name}...")

    image_original = imread(tif_path)

    if image_original.ndim != 2:
        raise ValueError(
            f"Expected a 2D image, but {tif_path.name} has shape "
            f"{image_original.shape}."
        )

    if scale_factor != 1.0:
        image_rescaled = resize(
            image_original,
            (
                int(image_original.shape[0] * scale_factor),
                int(image_original.shape[1] * scale_factor),
            ),
            preserve_range=True,
            anti_aliasing=True,
        ).astype("float32")
    else:
        image_rescaled = image_original.astype("float32")

    gc.collect()

    if gpu_available:
        torch.cuda.empty_cache()

    print("Running Cellpose segmentation...")

    masks, flows, styles = model.eval(
        [image_rescaled],
        channels=[0, 0],
        diameter=diameter_rescaled,
        do_3D=False,
        augment=False,
        batch_size=1,
    )

    mask_rescaled = masks[0].astype("uint16")

    print(
        f"Segmentation complete. "
        f"Number of segmented objects: {int(mask_rescaled.max())}"
    )

    rescaled_output_path = output_dir / f"{tif_path.stem}_mask_rescaled.tif"

    imwrite(rescaled_output_path, mask_rescaled)

    if scale_factor != 1.0:
        mask_original_size = resize(
            mask_rescaled,
            image_original.shape,
            order=0,
            preserve_range=True,
            anti_aliasing=False,
        ).astype("uint16")
    else:
        mask_original_size = mask_rescaled

    original_size_output_path = output_dir / f"{tif_path.stem}_mask_original_size.tif"

    imwrite(original_size_output_path, mask_original_size)

    print(f"Saved: {rescaled_output_path}")
    print(f"Saved: {original_size_output_path}")
