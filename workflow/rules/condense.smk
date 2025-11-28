rule condense:
    input:
        script = os.path.join(script_dir, "condense.py"),
        table = os.path.join("results", "joint_data", "{ab}.feather")
    output:
        os.path.join("results", "condensed_data", "{ab}.pkl")
    conda:
        os.path.join(env_dir, "predictor.yaml")
    threads: 8
    shell:
        "python3 {input.script} {input.table}"
