import nibabel as nib
import numpy as np
import torch
from monai.networks.nets import UNet  # Optional, use your own if preferred
from lib.utils import setup_logger

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)

device = snakemake.params.device
model = UNet(
    spatial_dims=3,
    in_channels=1,
    out_channels=19,  # for 18 labels + background
    channels=(32, 64, 128, 256),
    strides=(2, 2, 2),
    num_res_units=1,
)

checkpoint = torch.load(snakemake.params.model_weights, map_location=device)
model.load_state_dict(checkpoint["model_state_dict"])
model.to(device)
model.eval()

# Run inference
nii = nib.load(snakemake.input.nii)
img = nii.get_fdata()
img_tensor = torch.from_numpy(img).float().unsqueeze(0).unsqueeze(0).to(device)
# shape becomes: [1, 1, D, H, W]

with torch.no_grad():
    output = model(img_tensor)
    pred = output.argmax(dim=1).squeeze(0).cpu().numpy().astype(np.uint8)

# Save the output
output_img = nib.Nifti1Image(pred, affine=nii.affine)
nib.save(output_img, snakemake.output.nii)
