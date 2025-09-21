rule condense:
    input:
        script=os.path.join("scripts","condense.py"),
        table=os.path.join("results","joint_data", "{ab}.feather")
    output:
        os.path.join("results","condensed_data", "{ab}.pkl")
    conda:
        "predictor"
    threads: 8
    shell:
        "python3 {input.script} {input.table}"
