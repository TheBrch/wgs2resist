rule train:
    input:
        script = os.path.join(script_dir, "train.py"),
        src = os.path.join("results", "condensed_data", "{ab}.pkl"),
        lab = os.path.join("results", "binarized_data", "{ab}_lab.pkl"),
        sus = suscept
    output:
        expand(
            os.path.join("results", "models", "{{ab}}", "{model}.pkl"),
            model=["gaussian", "svm", "logistic"],
        )
    conda:
        os.path.join(env_dir, "predictor.yaml")
    threads: 4
    shell:
        "python3 {input.script} {input.src} {input.lab}"
