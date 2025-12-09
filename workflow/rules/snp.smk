rule consensus:
    input:
        aln=os.path.join(gene_aln_dir, "{gene}.aln.fas")
    output:
        consensus=os.path.join("results", "consensus_genes", "{gene}.fa")
    conda:
        os.path.join(env_dir, "emboss.yaml")
    params:
        custom_header=lambda wildcards, input: f">{os.path.basename(input.aln).split(".")[0]}"
    shell:
        "cons -sequence {input} -outseq /dev/stdout | sed '1c\\{params.custom_header}' > {output.consensus} "

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
        gff=os.path.join("results", "bakta", "reference", "reference.gff")
    threads: 8
    conda:
        os.path.join(env_dir, "bakta.yaml")
    params:
        gff_dir = subpath(output.gff, parent=True),
        db_dir = subpath(input.database, parent=True)
    shell:
        "bakta {input.fasta} --output {params.gff_dir} --db {params.db_dir} --force"
    

rule snippy:
    input:
        a=os.path.join("resources", "fasta_symlinks", "{id}.fasta"),
        ref=os.path.join("results", "bakta", "reference", "reference.gff")
    output:
        tab=os.path.join("results", "snippy", "{id}", "snps.tab")
    conda:
        os.path.join(env_dir, "snippy.yaml")
    params:
        soft=config["snippy"]["maxsoft"],
        frac=config["snippy"]["minfrac"],
        cov=config["snippy"]["mincov"],
        outdir=lambda wildcards, output: subpath(output.tab)
    shell:
        "snippy --outdir {output.aln} --maxsoft {params.soft} --minfrac {params.frac} --mincov {params.cov} --ref {input.ref} --ctgs {input.a} --force"

# rule core:
#     input:
#         expand(os.path.join("results", "snippy", "{sample}", "snps.tab"), sample=ids)
#     output:
#         core=os.path.join("results", "snippy-core", "core.tab")
#     run:
#         with open(output, "w") as o:
#             for f in input:
#                 for line in f:
#                     fields = line.strip().split('\t')
#                     chrom=fields[0]
#                     pos=fields[1]
#                     typ=fields[2]
#                     ref=fields[3]
#                     alt=fields[4]
    # filtering out non-synonymous variants