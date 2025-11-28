rule figs:
    input:
        script = os.path.join(script_dir, "figs.R"),
        src = expand(
            os.path.join("results", "models", "{{ab}}", "{model}.pkl"),
            model=["gaussian", "svm", "logistic"],
        )
    output:
        os.path.join("results", "models", "{ab}", "stats", "figs", "{ab}_pr.png")
    conda:
        os.path.join(env_dir, "R.yaml")
    shell:
        "Rscript {input.script} {wildcards.ab}"
