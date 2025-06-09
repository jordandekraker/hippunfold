# Template-based segmentation supports templates that have only a single hemisphere
# by flipping it


rule prep_segs_for_greedy:
    input:
        "{prefix}_dseg.nii.gz",
    params:
        labels=" ".join(str(label) for label in config["shape_inject"]["labels_reg"]),
        smoothing_stdev=config["shape_inject"]["label_smoothing_stdev"],
    output:
        temp(directory("{prefix}_dsegsplit")),
    group:
        "subj"
    conda:
        "../envs/c3d.yaml"
    shell:
        "mkdir -p {output} && "
        "c3d {input} -retain-labels {params.labels} -split -foreach -smooth {params.smoothing_stdev} -endfor -oo {output}/label_%02d.nii.gz"


def get_image_pairs(wildcards, input):
    """This rule requires snakemake 6.4.0, since it uses the new feature to execute if input files are not found"""

    import errno

    if not os.path.exists(input.subject_seg):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), input.subject_seg
        )
    if not os.path.exists(input.template_seg):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), input.template_seg
        )

    args = []
    # prep_segs_for_greedy creates label_{i} images for each entry in labels_reg,
    # but the numbering will be from 1 to N (not the numbers in the list)
    for label in range(1, len(config["shape_inject"]["labels_reg"]) + 1):
        subject_label = f"{input.subject_seg}/label_{label:02d}.nii.gz"
        template_label = f"{input.template_seg}/label_{label:02d}.nii.gz"

        if not os.path.exists(subject_label):
            print(f"Warning: {subject_label} does not exist, not using in registration")
            continue
        if not os.path.exists(template_label):
            print(
                f"Warning: {template_label} does not exist, not using in registration"
            )
            continue

        args.append("-i")
        args.append(subject_label)  # subject is fixed
        args.append(template_label)  # template is moving
    return " ".join(args)


def copy_or_flip(wildcards, file_to_process):
    if wildcards.hemi == "R":
        cmd = f"cp {file_to_process}"
    else:
        cmd = f"c3d {file_to_process} -flip x -o"
    return cmd


rule import_template_dseg:
    input:
        template=Path(workflow.basedir)
        / "../resources/upenn_layers/tpl-upenn_desc-hipptissue_layers.nii.gz",
    params:
        copy_or_flip_cmd=lambda wildcards, input: copy_or_flip(
            wildcards, input.template
        ),
    output:
        template_seg=temp(
            bids(
                root=root,
                datatype="anat",
                space="template",
                **inputs.subj_wildcards,
                desc="hipptissue",
                hemi="{hemi}",
                suffix="dseg.nii.gz",
            )
        ),
    group:
        "subj"
    conda:
        "../envs/c3d.yaml"
    shell:
        "{params.copy_or_flip_cmd} {output.template_seg}"


rule resample_template_dseg_tissue_for_reg:
    input:
        template_seg=bids(
            root=root,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            desc="hipptissue",
            hemi="{hemi}",
            suffix="dseg.nii.gz",
        ),
    params:
        resample_cmd="-resample-mm {res}".format(
            res=config["resample_dseg_for_templatereg"]
        ),
        crop_cmd="-trim 5vox",  #leave 5 voxel padding
    output:
        template_seg=temp(
            bids(
                root=root,
                datatype="anat",
                space="template",
                **inputs.subj_wildcards,
                desc="hipptissueresampled",
                hemi="{hemi}",
                suffix="dseg.nii.gz",
            )
        ),
    conda:
        "../envs/c3d.yaml"
    group:
        "subj"
    shell:
        "c3d {input} -int 0 {params.resample_cmd} {params.crop_cmd} -o {output}"


rule template_shape_lamareg:
    # this will provide a coarse registration and hopefully align the outer contours before finer registration with greedy
    input:
        template_seg=bids(
            root=root,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            desc="hipptissueresampled",
            hemi="{hemi}",
            suffix="dseg.nii.gz",
        ),
        subject_seg=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            desc="nnunet",
            space="corobl",
            hemi="{hemi}",
            suffix="dseg.nii.gz",
        ),
    params:
        general_opts="-d 3 -m SSD",
        affine_opts="-moments 2 -det 1",
        labelsmask="1 2 3 4 5 6 7 8",
    output:
        affine=temp(
            bids(
                root=root,
                **inputs.subj_wildcards,
                suffix="xfm.mat",
                datatype="warps",
                desc="lamareg",
                from_="template",
                to="subject",
                space="corobl",
                type_="itk",
                hemi="{hemi}",
            )
        ),
        warp=temp(
            bids(
                root=root,
                **inputs.subj_wildcards,
                suffix="xfm.nii.gz",
                datatype="warps",
                desc="lamareg",
                from_="template",
                to="subject",
                space="corobl",
                hemi="{hemi}",
            )
        ),
        out=temp(
            bids(
                root=root,
                **inputs.subj_wildcards,
                suffix="dseg.nii.gz",
                datatype="anat",
                desc="lamareg",
                space="corobl",
                hemi="{hemi}",
            )
        ),
    group:
        "subj"
    # conda:
    #     "../envs/lamareg.yaml"
    log:
        bids_log("template_shape_lamareg", **inputs.subj_wildcards, hemi="{hemi}"),
    shell:
        "lamar coregister --fixed {input.subject_seg} --moving {input.template_seg} --rev-affine {output.affine} --rev-warp-file {output.warp} --output {output.out} &> {log}"


def get_inject_scaling_opt(wildcards):
    """sets the smoothness of the greedy template shape injection deformation"""

    gradient_sigma = 1.732 * float(config["inject_template_smoothing_factor"])
    warp_sigma = 0.7071 * float(config["inject_template_smoothing_factor"])

    return f"-s {gradient_sigma}vox {warp_sigma}vox"


rule template_shape_reg:
    input:
        template_seg=bids(
            root=root,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            desc="hipptissueresampled",
            hemi="{hemi}",
            suffix="dsegsplit",
        ),
        subject_seg=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dsegsplit",
            desc="nnunet",
            space="corobl",
            hemi="{hemi}",
        ),
        affine=bids(
            root=root,
            **inputs.subj_wildcards,
            suffix="xfm.mat",
            datatype="warps",
            desc="lamareg",
            from_="template",
            to="subject",
            space="corobl",
            type_="ras",
            hemi="{hemi}",
        ),
        warp=bids(
            root=root,
            **inputs.subj_wildcards,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="lamareg",
            from_="template",
            to="subject",
            space="corobl",
            hemi="{hemi}",
        ),
    params:
        general_opts="-d 3 -m SSD",
        affine_opts="-moments 2 -det 1",
        greedy_opts=get_inject_scaling_opt,
        img_pairs=get_image_pairs,
    output:
        warp=temp(
            bids(
                root=root,
                **inputs.subj_wildcards,
                suffix="xfm.nii.gz",
                datatype="warps",
                desc="greedy",
                from_="template",
                to="subject",
                space="corobl",
                hemi="{hemi}",
            )
        ),
    group:
        "subj"
    conda:
        "../envs/greedy.yaml"
    threads: 8
    log:
        bids_log("template_shape_reg", **inputs.subj_wildcards, hemi="{hemi}"),
    shell:
        "greedy -threads {threads} {params.general_opts} {params.greedy_opts} {params.img_pairs} -ia {input.affine} -it {input.warp} -o {output.warp} &>> {log}"


rule template_shape_inject:
    input:
        template_seg=bids(
            root=root,
            datatype="anat",
            space="template",
            **inputs.subj_wildcards,
            desc="hipptissue",
            hemi="{hemi}",
            suffix="dseg.nii.gz",
        ),
        ref=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="nnunet",
            space="corobl",
            hemi="{hemi}",
        ),
        warp=bids(
            root=root,
            **inputs.subj_wildcards,
            suffix="xfm.nii.gz",
            datatype="warps",
            desc="greedy",
            from_="template",
            to="subject",
            space="corobl",
            hemi="{hemi}",
        ),
    params:
        interp_opt="-ri LABEL 0.1mm",  # smoothing sigma = 100micron
    output:
        inject_seg=temp(
            bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                suffix="dseg.nii.gz",
                desc="inject",
                space="corobl",
                hemi="{hemi}",
            )
        ),
    log:
        bids_log(
            "template_shape_inject",
            **inputs.subj_wildcards,
            hemi="{hemi}",
        ),
    group:
        "subj"
    conda:
        "../envs/greedy.yaml"
    threads: 8
    shell:
        "greedy -d 3 -threads {threads} {params.interp_opt} -rf {input.ref} -rm {input.template_seg} {output.inject_seg} -r {input.warp} &> {log}"


rule reinsert_subject_labels:
    """ c3d command to:
                1) get the labels to retain
                2) reslice to injected seg
                3) limit to labels_overwrite (e.g. SRLM)
                4) set injected seg to zero where retained labels are
                5) and add this to retained labels"""
    input:
        inject_seg=bids(
            root=root,
            datatype="anat",
            **inputs.subj_wildcards,
            suffix="dseg.nii.gz",
            desc="inject",
            space="corobl",
            hemi="{hemi}",
        ),
        subject_seg=get_input_for_shape_inject,
    params:
        labels=" ".join(
            str(label) for label in config["shape_inject"]["labels_reinsert"]
        ),
        labels_overwrite=" ".join(
            str(label) for label in config["shape_inject"]["labels_overwrite"]
        ),
    output:
        postproc_seg=temp(
            bids(
                root=root,
                datatype="anat",
                **inputs.subj_wildcards,
                suffix="dseg.nii.gz",
                desc="postproc",
                space="corobl",
                hemi="{hemi}",
            )
        ),
    group:
        "subj"
    conda:
        "../envs/c3d.yaml"
    shell:
        "c3d {input.subject_seg} -retain-labels {params.labels} -popas LBL "
        " -int 0 {input.inject_seg} -as SEG -push LBL -reslice-identity -popas LBL_RESLICE "
        " -push SEG -retain-labels {params.labels_overwrite}  -binarize "
        " -push LBL_RESLICE -multiply -popas LBL_RESLICE "
        "-push LBL_RESLICE -threshold 0 0 1 0 -push SEG -multiply "
        "-push LBL_RESLICE -add -o {output.postproc_seg}"
