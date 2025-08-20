import numpy as np


def get_gm_labels(wildcards):
    lbl_list = " ".join(
        [str(lbl) for lbl in config["laplace_labels"][wildcards.label]["gm"]]
    )
    return lbl_list


def get_src_sink_labels(wildcards):
    lbl_list = " ".join(
        [
            str(lbl)
            for lbl in config["laplace_labels"][wildcards.label][wildcards.dir][
                wildcards.srcsink
            ]
        ]
    )
    return lbl_list


def get_nan_labels(wildcards):
    gm_labels = set(map(int, config["laplace_labels"][wildcards.label]["gm"]))
    all_labels = set(np.arange(1, 18))
    missing_labels = sorted(all_labels - gm_labels)
    lbl_list = " ".join([str(lbl) for lbl in missing_labels])
    return lbl_list


rule get_label_mask:
    input:
        seg=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="postproc",
            space="corobl",
            hemi="{hemi}",
        ),
    params:
        labels=get_gm_labels,
    output:
        mask=temp(
            bids(
                root=root,
                datatype="coords",
                suffix="mask.nii.gz",
                space="corobl",
                desc="GM",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    group:
        "subj"
    conda:
        "../envs/c3d.yaml"
    shell:
        "c3d {input} -background -1 -retain-labels {params} -binarize {output}"


rule get_src_sink_mask:
    input:
        seg=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="postproc",
            space="corobl",
            hemi="{hemi}",
        ),
    params:
        labels=get_src_sink_labels,
    output:
        mask=temp(
            bids(
                root=root,
                datatype="coords",
                suffix="mask.nii.gz",
                space="corobl",
                dir="{dir}",
                desc="{srcsink,src|sink}",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    conda:
        "../envs/c3d.yaml"
    group:
        "subj"
    shell:
        "c3d {input} -background -1 -retain-labels {params} -binarize {output}"


rule get_src_sink_sdt:
    """calculate signed distance transform (negative inside, positive outside)"""
    input:
        mask=bids(
            root=root,
            datatype="coords",
            suffix="mask.nii.gz",
            space="corobl",
            dir="{dir}",
            desc="{srcsink}",
            hemi="{hemi}",
            label="{label}",
            **inputs.subj_wildcards,
        ),
    output:
        sdt=temp(
            bids(
                root=root,
                datatype="coords",
                suffix="sdt.nii.gz",
                space="corobl",
                dir="{dir}",
                desc="{srcsink}",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    conda:
        "../envs/c3d.yaml"
    group:
        "subj"
    shell:
        "c3d {input} -sdt -o {output}"


rule get_nan_mask:
    input:
        seg=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="postproc",
            space="corobl",
            hemi="{hemi}",
        ),
    params:
        labels=get_nan_labels,
    output:
        mask=temp(
            bids(
                root=root,
                datatype="coords",
                suffix="mask.nii.gz",
                space="corobl",
                dir="{dir}",
                desc="nan",
                hemi="{hemi}",
                label="{label}",
                **inputs.subj_wildcards,
            )
        ),
    conda:
        "../envs/c3d.yaml"
    group:
        "subj"
    shell:
        "c3d {input} -background -1 -retain-labels {params} -binarize {output}"


rule smooth_synthlayer:
    input:
        dseg_tissue=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="postproc",
            space="corobl",
            hemi="{hemi}",
        ),
    params:
        gm_labels=lambda wildcards: " ".join(
            [str(lbl) for lbl in config["laplace_labels"][wildcards.label]["gm"]]
        ),
        src_labels=lambda wildcards: " ".join(
            [
                str(lbl)
                for lbl in config["laplace_labels"][wildcards.label][wildcards.dir][
                    "src"
                ]
            ]
        ),
        sink_labels=lambda wildcards: " ".join(
            [
                str(lbl)
                for lbl in config["laplace_labels"][wildcards.label][wildcards.dir][
                    "sink"
                ]
            ]
        ),
        sigma=1.0,  # in voxels
    output:
        equidist=temp(
            bids(
                root=root,
                datatype="coords",
                dir="{dir,IO}",
                label="{label}",
                suffix="coords.nii.gz",
                desc="equidist",
                space="corobl",
                hemi="{hemi}",
                **inputs.subj_wildcards,
            )
        ),
    conda:
        "../envs/pyunfold.yaml"
    log:
        bids_log(
            "smooth_synthlayer",
            **inputs.subj_wildcards,
            dir="{dir, IO}",
            label="{label}",
            hemi="{hemi}",
        ),
    group:
        "subj"
    script:
        "../scripts/smooth_synthlayer.py"
