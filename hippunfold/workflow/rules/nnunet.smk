rule run_inference:
    """ This rule uses either GPU or CPU .
    It also runs in an isolated folder (shadow), with symlinks to inputs in that folder, copying over outputs once complete, so temp files are not retained"""
    input:
        nii=bids(
            root=root,
            datatype="anat",
            space="corobl",
            hemi="{hemi,L|R}",
            suffix="preproc.nii.gz",
            **inputs.subj_wildcards,
        ),
    params:
        model_weights=workflow.basedir + "/../resources/models/model_epoch2k.ckpt",
        device="cuda" if config["use_gpu"] else "cpu",
    output:
        nii=temp(
            bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                suffix="dseg.nii.gz",
                desc="nnunet",
                space="corobl",
                hemi="{hemi}",
            )
        ),
    log:
        bids_log(
            "run_inference",
            **inputs.subj_wildcards,
            hemi="{hemi}",
        ),
    threads: 16
    resources:
        gpus=1 if config["use_gpu"] else 0,
        mem_mb=16000,
        time=30 if config["use_gpu"] else 60,
    group:
        "subj"
    conda:
        "../envs/torch.yaml"
    script:
        "../scripts/torch_inference.py"


rule qc_nnunet_dice:
    input:
        res_mask=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="nnunet",
            space="corobl",
            hemi="{hemi}",
        ),
    params:
        hipp_lbls=[1, 2, 3, 4, 5, 6, 7, 8],
        ref=lambda wildcards: (
            Path(workflow.basedir)
            / "../resources/CITI168-slim/Mask_200umCoronalOblique_hemi-{hemi}.nii.gz".format(
                **wildcards
            )
        ),
    output:
        dice=report(
            bids(
                root=root,
                datatype="qc",
                suffix="dice.tsv",
                desc="unetf3d",
                hemi="{hemi}",
                **inputs.subj_wildcards,
            ),
            caption="../report/nnunet_qc.rst",
            category="Segmentation QC",
        ),
    group:
        "subj"
    conda:
        "../envs/pyunfold.yaml"
    script:
        "../scripts/dice.py"
