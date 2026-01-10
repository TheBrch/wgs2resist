rule unitig_input:
    output:
        os.path.join("results", "unitig-caller", "fastafiles.txt")
    run:
        with open(output[0], "w") as f:
            f.writelines(f"{p}\n" for p in sorted(fastapaths.values()))

rule call_unitigs:
    input:
        refs=os.path.join("results", "unitig-caller", "fastafiles.txt")
    output:
        os.path.join("results", "unitig-caller", "unitigs.rtab")
    conda:
        os.path.join(env_dir, "unitig-caller.yaml")
    threads:
        workflow.cores
    params:
        outdir=lambda wildcards, output: os.path.join(os.path.dirname(output[0]), os.path.basename(output[0]).split('.')[0])
    shell:
        '''
        unitig-caller --call --refs {input.refs} --rtab --write-graph --threads {threads} --out {params.outdir}
        '''
    
rule filter_unitigs:
    input:
        script = os.path.join(script_dir, "filter_unitigs.py"),
        uni = os.path.join("results", "unitig-caller", "unitigs.rtab"),
        meta = os.path.abspath(metadata)
    output:
        os.path.join("results", "unitig-caller", "filtered_unitigs.rtab")
    conda:
        os.path.join(env_dir, "predictor.yaml")
    shell:
        'python3 {input.script} {input.uni} {input.meta} > {output}'