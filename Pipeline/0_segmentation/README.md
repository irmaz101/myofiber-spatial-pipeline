# Segmentation

Myofibres and nuclei are segmented separately with Cellpose-SAM 4.0.4. The
scripts support patch annotation, custom-model training and full-image
segmentation with either the default or a trained model.

## Environment

```bash
conda create -n myofiber-segmentation python=3.10
conda activate myofiber-segmentation
```

Install the PyTorch build for your hardware using the
[PyTorch instructions](https://pytorch.org/get-started/locally/), then install:

```bash
python -m pip install -r Pipeline/0_segmentation/requirements.txt
```

## Images and settings

Inputs are two-dimensional TIFF images. The study used the combined
ATP1A1/E-cadherin/CD45 membrane signal for myofibres and DAPI for nuclei.
Use separate image, label and model directories for the two tasks.

| Setting | Value |
| --- | --- |
| Base model | Cellpose-SAM |
| Training patches | 12 per segmentation task |
| Patch size | 1800 × 1800 pixels |
| Training epochs | 800 |
| Full-image scale factor | 0.5 |
| Background label | 0 |

## Run the steps

Run each command from the repository root after setting the paths in the
corresponding script.

### 1. Extract patches

```bash
python Pipeline/0_segmentation/0_patches.py
```

Images are read from `raw_images/` and patches are saved in
`training_patches/`. The default `patches_per_image = 12` generates 12 patches
for every input image. Adjust this value or the image selection when using
multiple images to obtain the study total of 12 patches per task. Patch
coordinates are included in filenames; the sampling seed is 0.

### 2. Annotate patches

```bash
python Pipeline/0_segmentation/1_annotation_of_patches.py
```

Set `patch_dir`, `label_dir` and `sample_pattern`. With `annotate = True`,
patches open in Napari and masks are saved as `*_label.png` in
`training_labels/`. Assign a different positive integer label to each object
and leave background as `0`. With `annotate = False`, the patch locations are
shown on the original images instead.

### 3. Train a model

```bash
python Pipeline/0_segmentation/2_train_cellpose_model.py
```

Training matches each `*_label.png` mask to its TIFF patch and fine-tunes
Cellpose-SAM for 800 epochs. Set `patch_dir`, `label_dir`, `model_dir` and
`model_name` before running. The study trained separate fibre and nuclear
models. Training requests a GPU; full-image segmentation checks CUDA
availability automatically.

### 4. Segment full images

```bash
python Pipeline/0_segmentation/3_segment_objects.py
```

Set `input_dir` and `output_dir`. Leave `model_path = None` for the default
model, or provide the path to trained weights. Images are downscaled by
`scale_factor = 0.5`; a supplied `diameter_orig` is scaled accordingly.

Two masks are saved per image:

| Output | Dimensions |
| --- | --- |
| `<sample>_mask_rescaled.tif` | Image dimensions used for segmentation |
| `<sample>_mask_original_size.tif` | Original image dimensions, restored by nearest-neighbour resizing |

Use the original-size mask for transcript assignment. Inspect it for missing,
merged or fragmented objects, and confirm that it matches the transcript
coordinate system. Before preprocessing, copy or rename it to match the
transcript basename, for example
`cellpose_output/sample_01_mask_original_size.tif` to `masks_mf/sample_01.tif`.
