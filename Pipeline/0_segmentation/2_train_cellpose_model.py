"""Fine-tune Cellpose-SAM on matching image patches and manual masks."""

from pathlib import Path

from cellpose import io, models, train

# Settings

patch_dir = Path("training_patches")
label_dir = Path("training_labels")
model_dir = Path("models")

n_epochs = 800
model_name = f"cellpose_{n_epochs}epochs"

model_dir.mkdir(parents=True, exist_ok=True)


# Load training data

patches = []
labels = []

print("Loading training patches and labels...")

for label_path in sorted(label_dir.glob("*_label.png")):

    patch_name = label_path.name.replace("_label.png", ".tif")
    patch_path = patch_dir / patch_name

    if not patch_path.exists():
        raise FileNotFoundError(
            f"No matching image patch found for {label_path.name}: " f"{patch_path}"
        )

    patch = io.imread(patch_path)
    label = io.imread(label_path)

    patches.append(patch)
    labels.append(label)

    print(f"Matched {patch_path.name} with {label_path.name}")


if not patches:
    raise FileNotFoundError(
        f"No labelled training patches were found in {label_dir.resolve()}."
    )

print(f"Loaded {len(patches)} image-label pairs.")


# Initialize and train Cellpose

model = models.CellposeModel(gpu=True)

trained_model_path = train.train_seg(
    model.net,
    train_data=patches,
    train_labels=labels,
    n_epochs=n_epochs,
    model_name=model_name,
    save_path=model_dir,
    min_train_masks=1,
)

print(f"Training complete. Model saved to: {trained_model_path}")
