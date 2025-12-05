rule figs:
    input:
        script = os.path.join(script_dir, "figs.R"),
        src = expand(
            os.path.join("results", "models", "{{ab}}", "{model}.pkl"),
            model=models,
        )
    output:
        os.path.join("results", "models", "{ab}", "stats", "figs", "{ab}_pr.png")
    conda:
        os.path.join(env_dir, "R.yaml")
    shell:
        "Rscript {input.script} {wildcards.ab}"

rule summary_fig:
    input:
        script = os.path.join(script_dir, "fig_summary.R"),
        files = lambda x: expand(
            os.path.join("results", "models", "{ab}", "stats", "{model}_crossval_results.tsv"),
            ab=get_antibiotics(),
            model=models
        )
    output:
        os.path.join("results", "models", "combined_metrics_summary.png")
    conda:
        os.path.join(env_dir, "R.yaml")
    shell:
        "Rscript {input.script}"