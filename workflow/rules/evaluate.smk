rule figs:
    input:
        script = os.path.join(script_dir, "figs.R"),
        crossval = expand(
            os.path.join("results", "models", "{{ab}}", "stats", "{model}_crossval_results.tsv"),
            model=models
        ),
        svn_annot = os.path.join("results", "ns-snippy-core", "snv-effects.tsv"),
        gpa_annot = os.path.join("results", "bakta", "pan_genome_reference", "pan_genome_reference.tsv")
    output:
        os.path.join("results", "models", "{ab}", "stats", "figs", "{ab}_pr.png")
    conda:
        os.path.join(env_dir, "R.yaml")
    threads: workflow.cores/2
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