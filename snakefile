import glob, os

if os.path.exists('config.yaml'):
    configfile: 'config.yaml'

source=config['source_table']
suscept=config['suscept_table']

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
        lambda x: expand("models/{ab}/{model}.pkl", ab=get_antibiotics(x), model=["gaussian","svm","logistic","xgboost"]),

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
        "binarized_data/{ab}.pkl",
        "binarized_data/{ab}_lab.pkl"
    conda:
        "predictor"
    shell:
        "python3 {input.script} {input.table}"

rule condense:
    input:
        script="condense.py",
        table="binarized_data/{ab}.pkl"
    output:
        "condensed_data/{ab}.pkl"
    conda:
        "predictor"
    shell:
        "python3 {input.script} {input.table}"

rule train:
    input:
        script="train.py",
        src="condensed_data/{ab}.pkl",
        lab="binarized_data/{ab}_lab.pkl",
        sus=suscept
    output:
        expand("models/{{ab}}/{model}.pkl", model=["gaussian", "svm", "logistic", "xgboost"])
    conda:
        "predictor"
    shell:
        "python3 {input.script} {input.src} {input.lab}"