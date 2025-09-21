rule figs:
    input:
        script=os.path.join("scripts","figs.R"),
        src=expand(os.path.join("results","models", "{{ab}}", "{model}.pkl"), model=["gaussian", "svm", "logistic"]),
    output:
        os.path.join("results","models", "{ab}", "stats", "figs", "{ab}_pr.png")
    conda:
        "R"
    shell:
        "Rscript {input.script} {wildcards.ab}"