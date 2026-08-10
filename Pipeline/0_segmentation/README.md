# Segmentation

This directory contains the image-segmentation workflow used to generate labelled myofibre and nuclear masks from Xenium morphology images using Cellpose-SAM (Cellpose v4.0.4).

Myofibres and nuclei were segmented separately because skeletal muscle fibres are large multinucleated structures and cannot be represented adequately using nucleus-based cell segmentation.

## Files

| File                       | Purpose                                                                                     |
| -------------------------- | ------------------------------------------------------------------------------------------- |
| `0_patches.py`               | Extracts random image patches for manual annotation and model training.                     |
| `1_annotation_of_patches.py` | Opens the extracted patches in Napari for manual object annotation.                         |
| `2_train_cellpose_model.py`  | Trains a custom Cellpose-SAM model using manually annotated patches.                        |
| `3_segment_objects.py`       | Applies either the default Cellpose-SAM model or a custom trained model to complete images. |

The same scripts can be used for myofibre and nuclear segmentation by changing the input images, model path, and expected object diameter.

## Installation

Create a dedicated environment from the repository root:

```bash
conda create -n myofiber-segmentation python=3.10
conda activate myofiber-segmentation
```

Install a PyTorch build appropriate for the available CPU or CUDA platform
following the [official PyTorch instructions](https://pytorch.org/get-started/locally/),
then install the remaining dependencies:

```bash
python -m pip install -r Pipeline/0_segmentation/requirements.txt
```


## Input images

The scripts expect two-dimensional TIFF images.

The following Xenium morphology channels were used:

* Myofibres: combined ATP1A1/E-cadherin/CD45 membrane signal
* Nuclei: DAPI signal

Input paths and analysis parameters are defined in the **Settings** section near the beginning of each script.

## Workflow

### 1. Extract training patches

```bash
python Pipeline/0_segmentation/0_patches.py
```

For each segmentation task, 12 randomly selected `1800 × 1800` pixel patches
were extracted for manual annotation. Cellpose may internally sample smaller
tiles from these annotated regions during training.

### 2. Annotate training patches

```bash
python Pipeline/0_segmentation/1_annotation_of_patches.py
```

The extracted patches are opened in Napari for manual annotation.

Each object should be assigned a unique positive integer label, while background pixels should remain labelled as `0`. The annotated masks must correspond to their original image patches.

### 3. Train a custom Cellpose-SAM model

```bash
python Pipeline/0_segmentation/2_train_cellpose_model.py
```

Separate custom models were trained for myofibres and nuclei. Both models were initialized from the pretrained Cellpose-SAM model and trained for 800 epochs.

Before training, configure:

* the directory containing training images;
* the directory containing the corresponding masks;
* the output directory;
* the model name;
* the number of training epochs.

### 4. Segment complete images

```bash
python Pipeline/0_segmentation/3_segment_objects.py
```

To use the default pretrained Cellpose-SAM model, set:

```python
model_path = None
```

To use a custom trained model, provide the path to the model weights:

```python
model_path = Path("models/custom_model")
```

## Image resizing

Complete images are downscaled before segmentation using:

```python
scale_factor = 0.5
```

The object diameter is adjusted automatically according to the scale factor.

After segmentation, the masks are resized to the original image dimensions using nearest-neighbour interpolation. This preserves the integer object labels.

## Outputs

For each input image, the script saves:

```text
<sample>_mask_rescaled.tif
<sample>_mask_original_size.tif
```

The rescaled mask corresponds to the image dimensions used during Cellpose segmentation.

The original-size mask has the same dimensions as the input image and should be used for transcript assignment and downstream analysis.

## Mask requirements

A valid output mask must:

* have the same dimensions as the corresponding transcript-coordinate image;
* contain non-negative integer labels;
* use `0` for background;
* assign a unique positive integer label to each segmented object.

The final masks should be inspected visually before transcript assignment to identify missing, merged, fragmented, or incorrectly detected objects.

## Main study settings

| Parameter          | Setting                          |
| ------------------ | -------------------------------- |
| Cellpose version   | 4.0.4                            |
| Base model         | Cellpose-SAM                     |
| Training patches   | 12 per segmentation task         |
| Patch dimensions   | `1800 × 1800` pixels             |
| Training epochs    | 800                              |
| Image scale factor | 0.5                              |
| Background label   | 0                                |
| GPU use            | Automatic when CUDA is available |