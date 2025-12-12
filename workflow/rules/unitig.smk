rule unitig_input:
    output:
        os.path.join("resources", "fastafiles.txt")
    run:
        with open(output[0], "w") as f:
            f.writelines(f"{p}\n" for p in sorted(fastapaths.values()))

rule call_unitigs:
    input:
        refs=os.path.join("resources", "fastafiles.txt")
    output:
        os.path.join("results", "unitig-caller", "dummy.dummy")
    conda:
        os.path.join(env_dir, "unitig-caller.yaml")
    threads:
        workflow.cores
    params:
        outdir=subpath(output[0], parent=True)
    shell:
        '''
        unitig-caller --call --refs {input.refs} --rtab --write-graph --threads {threads} --out {params.outdir} ;
        touch {output}
        '''
    