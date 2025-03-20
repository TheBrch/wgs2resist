if os.path.exists('config.yaml'):
    configfile: 'config.yaml'

source=config['source_table']
suscept=config['suscept_table']

rule all:
    input:
        "gaussian.pkl",
        "svm.pkl",
        "xgboost.pkl",
        "logistic.pkl"

rule train:
    input:
        script="train.py",
        src="training_matrix.tsv",
        sus=suscept,
    output:
        "gaussian.pkl",
        "svm.pkl",
        "xgboost.pkl",
        "logistic.pkl"
    conda:
        "predictor"
    shell:
        "python3 {input.script}"

rule reformat:
    input:
        script="reformat.R"
    output:
        "training_matrix.tsv"
    conda:
        "R"
    shell:
        "Rscript {input.script}"
