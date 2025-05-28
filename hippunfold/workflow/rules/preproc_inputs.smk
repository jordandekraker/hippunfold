# populate the HIPPUNFOLD_CACHE_DIR folder as needed
from lib import utils as utils

download_dir = utils.get_download_dir()


def get_inputs():  # TODO: swap this in
    if config["modality"] == "multires":
        return inputs["T1w"] + inputs["T2w"] + inputs["FLAIR"]
    else:
        return inputs[config["modality"]]


rule import_any_modality:
    input:
        inputs[config["modality"]].path,
    output:
        bids(
            root=root,
            datatype="anat",
            suffix=config["modality"] + ".nii.gz",
            **inputs[config["modality"]].wildcards,
        ),
    group:
        "subj"
    shell:
        "cp {input} {output}"


rule lamareg_to_template:
    input:
        img=bids(
            root=root,
            datatype="anat",
            suffix=config["modality"] + ".nii.gz",
            **inputs[config["modality"]].wildcards,
        ),
        template_img=Path(workflow.basedir)
        / "../resources/CITI168-slim/T1w_space-corobl_1mm.nii.gz",
        template_seg=Path(workflow.basedir)
        / "../resources/CITI168-slim/synthseg_space-corobl.nii.gz",
    output:
        affine=bids(
            root=root,
            datatype="warps",
            suffix="xfm.mat",
            from_=config["modality"],
            to="corobl",
            type_="itk",
            **inputs[config["modality"]].wildcards,
        ),
        invaffine=bids(
            root=root,
            datatype="warps",
            suffix="xfm.mat",
            from_="corobl",
            to=config["modality"],
            type_="itk",
            **inputs[config["modality"]].wildcards,
        ),
        warp=bids(
            root=root,
            datatype="warps",
            suffix="xfm.nii.gz",
            from_=config["modality"],
            to="corobl",
            type_="itk",
            **inputs[config["modality"]].wildcards,
        ),
        invwarp=bids(
            root=root,
            datatype="warps",
            suffix="xfm.nii.gz",
            from_="corobl",
            to=config["modality"],
            type_="itk",
            **inputs[config["modality"]].wildcards,
        ),
        out=bids(
            root=root,
            datatype="anat",
            **inputs.wildcards,
            suffix=config["modality"],
            space="template",
        ),
    shadow:
        "minimal"
    # conda:
    #     conda_env("lamar")
    group:
        "subj"
    log:
        bids_log(
            "lamareg_to_template",
            **inputs.subj_wildcards,
        ),
    shell:
        "lamar register --fixed {input.template_img} --fixed-parc {input.template_seg} --moving {input.img} --affine {output.affine} --inverse-affine {output.invaffine} --warpfield {output.warp} --inverse-warpfield {output.invwarp} --output {output.out} &> {log}"


rule apply_transforms:
    input:
        img=bids(
            root=root,
            datatype="anat",
            suffix=config["modality"] + ".nii.gz",
            **inputs[config["modality"]].wildcards,
        ),
        template_img=Path(workflow.basedir)
        / "../resources/CITI168-slim/T1w_space-corobl_1mm.nii.gz",
        affine=bids(
            root=root,
            datatype="warps",
            suffix="xfm.mat",
            from_=config["modality"],
            to="corobl",
            type_="itk",
            **inputs[config["modality"]].wildcards,
        ),
        warp=bids(
            root=root,
            datatype="warps",
            suffix="xfm.nii.gz",
            from_=config["modality"],
            to="corobl",
            type_="itk",
            **inputs[config["modality"]].wildcards,
        ),
    output:
        img=bids(
            root=root,
            datatype="anat",
            space="corobl",
            hemi="{hemi,L|R}",
            suffix=config["modality"] + ".nii.gz",
            **inputs[config["modality"]].wildcards,
        ),
    group:
        "subj"
    log:
        bids_log(
            "apply_transforms",
            **inputs.subj_wildcards,
            hemi="{hemi}",
        ),
    shell:
        "lamar apply-warp --affine {input.affine} --warp {input.warp} --moving {input.img} --reference {input.template_img} --output {output.img} &> {log}"


# TODO: refine registrations using the raw images after the above initializations


rule template_xfm_itk2ras:
    input:
        "{prefix}_type-itk_xfm.mat",
    output:
        "{prefix}_type-ras_xfm.mat",
    conda:
        conda_env("c3d")
    group:
        "subj"
    shell:
        "c3d_affine_tool -itk {input} -o {output}"


rule superres_inputs:
    input:
        bids(
            root=root,
            datatype="anat",
            space="corobl",
            hemi="{hemi,L|R}",
            suffix=config["modality"] + ".nii.gz",
            **inputs[config["modality"]].wildcards,
        ),
    output:
        bids(
            root=root,
            datatype="anat",
            space="corobl",
            hemi="{hemi,L|R}",
            suffix="preproc.nii.gz",
            **inputs.subj_wildcards,
        ),
    conda:
        conda_env("c3d")
    group:
        "subj"
    shell:
        """
        if [ $(echo {input} | wc -w) -eq 1 ]; then
            cp {input} {output}
        else
            c3d {input} -add -o {output}
        fi
        """
