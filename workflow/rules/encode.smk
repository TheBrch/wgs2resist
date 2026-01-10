checkpoint reformat:
    input:
        script = os.path.join(script_dir, "reformat.R"),
        sus = suscept,
        tab = os.path.join("results", "ns-snippy-core", "non-synonymous-core.tsv"),
        unitigs = os.path.join("results", "unitig-caller", "filtered_unitigs.rtab")
    output:
        directory(os.path.join("results", "training_data"))
    conda:
        os.path.join(env_dir, "R.yaml")
    threads:
        workflow.cores/2
    shell:
        "Rscript {input.script} {input.tab} {input.unitigs}"


rule binarize:
    input:
        script = os.path.join(script_dir, "binarize.py"),
        table = os.path.join("results", "training_data", "{ab}.tsv")
    output:
        os.path.join("results", "binarized_data", "{ab}.feather"),
        os.path.join("results", "binarized_data", "{ab}_lab.pkl"),
    conda:
        os.path.join(env_dir, "predictor.yaml")
    shell:
        "python3 {input.script} {input.table}"


rule join:
    input:
        script = os.path.join(script_dir, "snp_gpa_join.R"),
        table = os.path.join("results", "binarized_data", "{ab}.feather"),
        gpa = gpam
    output:
        os.path.join("results", "joint_data", "{ab}.feather")
    conda:
        os.path.join(env_dir, "R.yaml")
    shell:
        "Rscript {input.script} {input.table}"
