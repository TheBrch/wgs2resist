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
        protected(os.path.join("resources", "db", "bakta.db"))
    conda:
        os.path.join(env_dir, "bakta.yaml")
    threads: 2
    params:
        output_dir = subpath(output[0], ancestor=2)
    shell:
        "bakta_db download --output {params.output_dir}"

rule snv_annotate:
    input:
        fasta=os.path.join("results", "reference.fa"),
        database=os.path.join("resources", "db", "bakta.db")
    output:
        gbff=os.path.join("results", "bakta", "reference", "reference.gbff"),
        tsv=os.path.join("results", "bakta", "reference", "reference.tsv")
    threads: workflow.cores
    conda:
        os.path.join(env_dir, "bakta.yaml")
    params:
        out_dir = subpath(output.gbff, parent=True),
        db_dir = subpath(input.database, parent=True)
    shell:
        "bakta {input.fasta} --output {params.out_dir} --db {params.db_dir} --threads {threads} --force"

rule gpa_annotate:
    input:
        fasta=pan_ref,
        database=os.path.join("resources", "db", "bakta.db")
    output:
        tsv = os.path.join("results", "bakta", "pan_genome_reference", "pan_genome_reference.tsv")
    threads: workflow.cores
    conda:
        os.path.join(env_dir, "bakta.yaml")
    params:
        out_dir = subpath(output.tsv, parent=True),
        db_dir = subpath(input.database, parent=True)
    shell:
        "bakta {input.fasta} --output {params.out_dir} --db {params.db_dir} --threads {threads} --force"

rule snippy:
    input:
        a=lambda wildcards: fastapaths[wildcards.id],
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
        os.path.join("results", "ns-snippy-core", "dn_ds.tsv"),
        os.path.join("results", "ns-snippy-core", "snv-effects.tsv")
    conda:
        os.path.join(env_dir, "custom_py.yaml")
    params:
        outdir=subpath(output[0], parent = True)
    threads: 4
    shell:
        '''
        python3 {input.script} -o {params.outdir} -s {input.samples}
        '''