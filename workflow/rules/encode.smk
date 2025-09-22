checkpoint reformat:
    input:
        script = os.path.join("scripts", "reformat.R"),
        sus = suscept
    output:
        directory(os.path.join("results", "training_data"))
    conda:
        "R"
    shell:
        "Rscript {input.script}"


rule binarize:
    input:
        script = os.path.join("scripts", "binarize.py"),
        table = os.path.join("results", "training_data", "{ab}.tsv")
    output:
        os.path.join("results", "binarized_data", "{ab}.feather"),
        os.path.join("results", "binarized_data", "{ab}_lab.pkl"),
    conda:
        "predictor"
    shell:
        "python3 {input.script} {input.table}"


rule join:
    input:
        script = os.path.join("scripts", "snp_gpa_join.R"),
        table = os.path.join("results", "binarized_data", "{ab}.feather"),
        gpa = gpam
    output:
        os.path.join("results", "joint_data", "{ab}.feather")
    conda:
        "R"
    shell:
        "Rscript {input.script} {input.table}"
