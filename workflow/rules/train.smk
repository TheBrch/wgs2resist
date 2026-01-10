rule train:
    input:
        script = os.path.join(script_dir, "train.py"),
        src = os.path.join("results", "condensed_data", "{ab}{suffix}.pkl"),
        lab = os.path.join("results", "binarized_data", "{ab}_lab.pkl"),
        sus = suscept
    wildcard_constraints:
        suffix = "(_unitig)?"
    output:
        expand(
            os.path.join("results", "models", "{{ab}}{{suffix}}", "{model}.pkl"),
            model=models,
        ),
        crossval = expand(
            os.path.join("results", "models", "{{ab}}{{suffix}}", "stats", "{model}_crossval_results.tsv"),
            model=models
        )
    conda:
        os.path.join(env_dir, "predictor.yaml")
    threads: 4
    shell:
        "python3 {input.script} {input.src} {input.lab}"
