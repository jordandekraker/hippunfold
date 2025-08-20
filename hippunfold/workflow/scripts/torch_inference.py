import os
import nibabel as nib
import numpy as np
import torch
from lib.utils import setup_logger
from monai.networks.nets import UNet

log_file = snakemake.log[0] if getattr(snakemake, "log", None) else None
logger = setup_logger(log_file)

# Resolve device
device = torch.device(snakemake.params.device)

# ---- Build a plain model (must match training setup!) ----
def build_model():
    # Mirror your training config here
    return UNet(
        spatial_dims=3,
        in_channels=1,
        out_channels=19,              # 18 labels + background
        channels=(32, 64, 128, 256),
        strides=(2, 2, 2),
        num_res_units=1,
    )

def strip_prefixes(state_dict):
    """Remove common Lightning/module prefixes so a plain nn.Module can load them."""
    prefixes = ("model.", "net.", "unet.", "module.", "seg_net.", "network.")
    clean = {}
    for k, v in state_dict.items():
        kk = k
        for p in prefixes:
            if kk.startswith(p):
                kk = kk[len(p):]
        clean[kk] = v
    return clean

# Load checkpoint on CPU first (safer) and extract state_dict
ckpt = torch.load(snakemake.params.model_weights, map_location="cpu")
state_dict = ckpt.get("state_dict", ckpt)
state_dict = strip_prefixes(state_dict)

model = build_model()
missing, unexpected = model.load_state_dict(state_dict, strict=False)
if missing or unexpected:
    logger.warning(f"Loaded with missing keys: {missing} | unexpected keys: {unexpected}")

model.to(device).eval()

# Load input image
nii = nib.load(snakemake.input.nii)
img = nii.get_fdata().astype(np.float32)
img_tensor = torch.from_numpy(img).unsqueeze(0).unsqueeze(0).to(device)  # [1,1,D,H,W]

with torch.inference_mode():
    logits = model(img_tensor)
    pred = logits.argmax(dim=1).squeeze(0).to("cpu", dtype=torch.uint8).numpy()

# Save output (preserve affine/header)
out = nib.Nifti1Image(pred, affine=nii.affine, header=nii.header)
nib.save(out, snakemake.output.nii)
logger.info(f"Wrote {snakemake.output.nii}")
