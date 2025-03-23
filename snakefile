if os.path.exists('config.yaml'):
    configfile: 'config.yaml'

source=config['source_table']
suscept=config['suscept_table']

antibiotics = []
with open(suscept, "r") as f:
    for line in f.readlines()[1:]:
        antibiotics.append(
            line.strip().split("\t")[0].replace("/", "_")
        )

rule all:
    input:
        lambda x: expand("models/{ab}/{model}.pkl", ab=antibiotics, model=["gaussian","svm","logistic","xgboost"]),

rule reformat:
    input:
        script="reformat.R"
    output:
        lambda x: expand(f"training_data/{ab}.tsv", ab=antibiotics)
    conda:
        "R"
    shell:
        "Rscript {input.script}"

rule binarize:
    input:
        script="binarize.py"
        table="training_data/{ab}.tsv"
    output:
        "binarized_data/{ab}.pkl"
    conda:
        "predictor"
    shell:
        "python3 {input.script} {input.table}"

rule condense:
    input:
        script="condense.py"
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
        src="binarized_data/{ab}.pkl",
        sus=suscept
    output:
        lambda x: expand(f"{ab}/{model}.pkl", model=["gaussian","svm","logistic","xgboost"], ab=[wildcards.ab])
    conda:
        "predictor"
    shell:
        "python3 {input.script} {input.src}"