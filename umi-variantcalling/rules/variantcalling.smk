# rule vardict_single_mode:
#     input:
#         reference=config["reference"]['genome'],
#         regions=config["reference"]["region_file"],
#         bam=lambda wildcards: input_bam[wildcards.sample],
#     output:
#         vcf= "vardict/{sample}.vcf",
#     params:
#         extra="-k 0 -U",
#         bed_columns="-c 1 -S 2 -E 3 -g 4",  # Optional, default is -c 1 -S 2 -E 3 -g 4
#         allele_frequency_threshold=config["vardict"]["vaf"],  # Optional, default is 0.01
#     threads: 10
#     benchmark:
#         "benchmarks/vardict/{sample}.tsv"
#     resources:
#         runtime="4h",
#         threads=10
#     log:
#         "logs/varscan_{sample}_s_.log",
#     wrapper:
#         "v5.6.1/bio/vardict"

rule vardict:
    input:
        ref=config["reference"]['genome'],
        # you can have a list of samples here
        samples=lambda wildcards: input_bam[wildcards.sample],
        # optional BED file specifying chromosomal regions on which freebayes
        # should run, e.g. all regions that show coverage
        regions=config["reference"]["region_file"],
    output:
        "vardict/{sample}.vcf"
    params:
        extra="-c 1 -S 2 -E 3 -g 4",
        vaf=config["vardict"]["vaf"]
    benchmark:
        "benchmarks/vardict/{sample}.tsv"
    resources:
        runtime="3h",
        mem="50G",       
        threads=10
    conda:
        "../env/vardict.yaml"
    shell:
        r"""
        mkdir -p vardict
        vardict-java -k 0 -U -G {input.ref} -f {params.vaf} -N {wildcards.sample} -b {input.samples} {params.extra} {input.regions} | ${{CONDA_PREFIX}}/share/vardict-java*/bin/teststrandbias.R | ${{CONDA_PREFIX}}/share/vardict-java*/bin/var2vcf_valid.pl -N {wildcards.sample} -E -f {params.vaf} > {output}
        """

rule bcf_to_vcf:
    input:
        "vardict/{sample}.vcf"
    output:
        temp("vardict/{sample}.vcf.gz")
    params:
        extra=""  # optional parameters for bcftools view (except -o)
    benchmark:
        "benchmarks/bcf_to_vcf/{sample}.tsv"
    group:
        "vcfgroup"
    log:
        "logs/bcf_to_vcf/{sample}.log"
    resources:
        runtime="30m"
    wrapper:
        "v5.0.2/bio/bcftools/view"

rule bcftools_index:
    input:
        "vardict/{sample}.vcf.gz"
    output:
        temp("vardict/{sample}.vcf.gz.csi")
    params:
        extra=""  # optional parameters for bcftools index
    benchmark:
        "benchmarks/bcftools_index/{sample}.tsv"
    group:
        "vcfgroup"
    resources:
        runtime="30m"
    wrapper:
        "v5.0.2/bio/bcftools/index"
    #    "v1.0.0/bio/bcftools/index"

rule norm_vcf:
    input:
        vcf = "vardict/{sample}.vcf.gz",
        index = "vardict/{sample}.vcf.gz.csi",
        ref=config["reference"]['genome'],
    output:
        temp("vardict/{sample}.norm.vcf.gz")
    params:
        extra="", 
        uncompressed_bcf=True,
    benchmark:
        "benchmarks/norm_vcf/{sample}.tsv"
    group:
        "vcfgroup"
    log:
        "logs/norm_vcf/{sample}.log"
    resources:
        runtime="30m"
    wrapper:
        "v5.5.2/bio/bcftools/norm"


rule vcf_to_table:
    input:
        "vardict/{sample}.norm.vcf.gz"
    output:
        temp("vardict/{sample}.csv")
    params:
        #wd = config["general"]['work_dir'] + "variantcalling/",
        sd = config["general"]['snakedir'],
        rs = f"scripts/Vcf2Table.R",
    benchmark:
        "benchmarks/vcf_to_table/{sample}.tsv"
    group:
        "vcfgroup"
    resources:
        runtime="30m"
    conda:
        "../env/vcf2table.yaml",
    shell:
        r"""
        Rscript {params.sd}/{params.rs} {input} {output}
        """
        #Rscript {params.sd}/{params.rs} {params.wd}/{input} {params.wd}/{output}

rule table_to_anno:
    input:
        "vardict/{sample}.csv",
    output:
        temp("table/{sample}.csv")
    conda:
        "../env/Renv.yaml"
    threads:
        1
    log:
        "logs/table_to_anno/{sample}.log"
    benchmark:
        "benchmarks/table_to_anno/{sample}.tsv"
    group:
        "vcfgroup"
    resources:
        runtime="30m"
    params:
        sd = config["general"]['snakedir'],
        rs = f"scripts/SplitTable.R",
    shell:
        r"""
        Rscript {params.sd}/{params.rs} {input} {output}
        """
# rule table_to_anno:
#     input:
#         "vardict/{sample}.csv",
#     output:
#         "table/{sample}.csv"
#     threads:
#         1
#     resources:
#         time="00:30:00"
#     params:
#         wd = config["general"]['work_dir'] + "variantcalling/",
#         sd = config["general"]['snakedir'],
#         rs = f"scripts/SplitTable.R",
#     shell:
#         r"""
#         Rscript {params.sd}/{params.rs} {params.wd}/{input} {params.wd}/{output}
#         """
