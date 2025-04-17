import glob, os

if os.path.exists('config.yaml'):
    configfile: 'config.yaml'

source=config['source_table']
suscept=config['suscept_table']
gpam=config['gpam_table']

def get_antibiotics (x):
    reformat_dir = checkpoints.reformat.get().output[0]
    r_tsvs = glob.glob(reformat_dir+"/*.tsv")
    antibiotics = []
    for table in r_tsvs:
        filename = table.split('/')[-1]
        match = re.match(r"(.+)\.tsv", filename)
        if match:
            name_part = match.groups()[0]
            antibiotics.append(name_part)
    return antibiotics

rule all:
    input:
        lambda x: expand("models/{ab}/stats/figs/{ab}_prc.png", ab=get_antibiotics(x))#,
        # lambda x: expand("condensed_data/{ab}.pkl", ab=get_antibiotics(x))

checkpoint reformat:
    input:
        script="reformat.R",
        sus=suscept
    output:
        directory("training_data")
    conda:
        "R"
    shell:
        "Rscript {input.script}"

rule binarize:
    input:
        script="binarize.py",
        table="training_data/{ab}.tsv"
    output:
        "binarized_data/{ab}.feather",
        "binarized_data/{ab}_lab.pkl"
    conda:
        "predictor"
    shell:
        "python3 {input.script} {input.table}"

rule join:
    input:
        script="snp_gpa_join.R",
        table="binarized_data/{ab}.feather",
        gpa=gpam
    output:
        "joint_data/{ab}.feather"
    conda:
        "R"
    shell:
        "Rscript {input.script} {input.table}"

rule condense:
    input:
        script="condense.py",
        table="joint_data/{ab}.feather"
    output:
        "condensed_data/{ab}.pkl"
    conda:
        "predictor"
    threads: 8
    shell:
        "python3 {input.script} {input.table}"

rule train:
    input:
        script="train.py",
        src="condensed_data/{ab}.pkl",
        lab="binarized_data/{ab}_lab.pkl",
        sus=suscept
    output:
        expand("models/{{ab}}/{model}.pkl", model=["gaussian", "svm", "logistic"])
    conda:
        "predictor"
    threads: 4
    shell:
        "python3 {input.script} {input.src} {input.lab}"

rule figs:
    input:
        script="figs.R",
        src=expand("models/{{ab}}/{model}.pkl", model=["gaussian", "svm", "logistic"]),
    output:
        "models/{ab}/stats/figs/{ab}_prc.png"
    conda:
        "R"
    shell:
        "Rscript {input.script} {wildcards.ab}"