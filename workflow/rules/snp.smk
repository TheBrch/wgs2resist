rule consensus:
    input:
        aln=os.path.join(gene_aln_dir, "{gene}.aln.fas")
    output:
        consensus=os.path.join("results", "consensus_genes", "{gene}.fa")
    conda:
        os.path.join(env_dir, "emboss.yaml")
    params:
        custom_header=lambda wildcards, input: os.path.basename(input.aln).split('.')[0]
    shell:
        '''
        cons -plurality 0 -name {params.custom_header} -sequence {input} -outseq /dev/stdout | \
        degapseq -sequence /dev/stdin -outseq {output.consensus}
        '''

rule gene_concat:
    input:
        expand(os.path.join("results", "consensus_genes", "{gene}.fa"), gene=[os.path.basename(i).split('.')[0] for i in gene_files])
    output:
        os.path.join("results", "reference.fa")
    shell:
        "cat {input} > {output}"

rule baktadb:
    output:
        os.path.join("resources", "db", "bakta.db")
    conda:
        os.path.join(env_dir, "bakta.yaml")
    threads: 2
    params:
        output_dir = subpath(output[0], ancestor=2)
    shell:
        "bakta_db download --output {params.output_dir}"

rule bakta:
    input:
        fasta=os.path.join("results", "reference.fa"),
        database=os.path.join("resources", "db", "bakta.db")
    output:
        gbff=os.path.join("results", "bakta", "reference", "reference.gbff")
    threads: 8
    conda:
        os.path.join(env_dir, "bakta.yaml")
    params:
        gbff_dir = subpath(output.gbff, parent=True),
        db_dir = subpath(input.database, parent=True)
    shell:
        "bakta {input.fasta} --output {params.gbff_dir} --db {params.db_dir} --force"
    

rule snippy:
    input:
        a=lambda wildcards: (true_id := get_file(wildcards.id)) and os.path.join(assembly_dir, true_id, f"{true_id}_filtered.fasta"),
        ref=os.path.join("results", "bakta", "reference", "reference.gbff")
    output:
        tab=os.path.join("results", "snippy", "{id}", "snps.tab"),
        csv=os.path.join("results", "snippy", "{id}", "snps.csv")
    conda:
        os.path.join(env_dir, "snippy.yaml")
    params:
        soft=config["snippy"]["maxsoft"],
        frac=config["snippy"]["minfrac"],
        cov=config["snippy"]["mincov"],
        outdir=subpath(output.tab, parent=True)
    shell:
        "snippy --outdir {params.outdir} --maxsoft {params.soft} --minfrac {params.frac} --mincov {params.cov} --ref {input.ref} --ctgs {input.a} --force"

rule non_syn_core:
    input:
        samples=lambda wildcards: expand(os.path.join("results", "snippy", "{ids}", "snps.csv"), ids=ids),
        script=os.path.join(script_dir, "filtered_snp_core.py")
    output:
        os.path.join("results", "ns-snippy-core", "non-synonymous-core.tsv"),
        os.path.join("results", "ns-snippy-core", "dn_ds.tsv")
    conda:
        os.path.join(env_dir, "custom_py.yaml")
    params:
        outdir=subpath(output[0], parent = True)
    threads: 4
    shell:
        '''
        python3 {input.script} -o {params.outdir} -s {input.samples}
        '''