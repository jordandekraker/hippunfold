import numpy as np
import nibabel as nib
from scipy.ndimage import gaussian_filter
from lib.utils import setup_logger

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)
logger.info("imports")

lbl = nib.load(snakemake.input.dseg_tissue)
newimg = np.zeros_like(lbl.get_fdata(), dtype=np.float32)
logger.info("loaded img")

gm_labels = [int(x) for x in snakemake.params.gm_labels.split()]
sink_labels = [int(x) for x in snakemake.params.sink_labels.split()]
source_labels = [int(x) for x in snakemake.params.src_labels.split()]
logger.info("loaded label lists")

for i in gm_labels:
    newimg[lbl.get_fdata() == i] = (i - min(gm_labels) + 1) / (len(gm_labels)+1)
for i in sink_labels:
    newimg[lbl.get_fdata() == i] = 1
logger.info("initial IO guess")


newimg = gaussian_filter(newimg, snakemake.params.sigma)
logger.info("smooth")

outimg = np.zeros_like(newimg, dtype=np.float32)
for i in gm_labels:
    outimg[lbl.get_fdata() == i] = newimg[lbl.get_fdata() == i]
logger.info("reset boundaries")

newimg = nib.Nifti1Image(newimg, lbl.affine, lbl.header)
nib.save(newimg, snakemake.output.equidist)