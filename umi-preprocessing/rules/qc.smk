# rule get_fastqc_input:
#     input:
#         lookup(query="sample == '{sample}'", within=SAMPLES, cols="full_path"),
#     output:
#         read1=temp("reads/{sample}.R1.fastq"),
#         read2=temp("reads/{sample}.R2.fastq")
#     group:
#         "fastqc"    
#     conda:
#             "../envs/mapping.yaml"
#     params:
#         musage=config["picard"]["memoryusage"]
#     benchmark:
#         "benchmarks/get_fastqc_input/{sample}.tsv"
#     log:
#         "logs/reads/{sample}.samtofastq.log"
#     threads:
#         8
#     resources:
#         mem_mb=get_mem_40_10,
#         runtime="1h"
#     shell:
#         r"""
#         picard {params.musage} SamToFastq I={input} F={output.read1} SECOND_END_FASTQ={output.read2} &> {log}
#         """

# rule fastqc:
#     input:
#         "reads/{sample}.{read}.fastq"
#     output:
#         html=temp("qc/fastqc/{sample}.{read}.html"),
#         zip=temp("qc/fastqc/{sample}.{read}_fastqc.zip") # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
#     group:
#         "fastqc" 
#     params: ""
#     benchmark:
#         "benchmarks/fastqc/{sample}.{read}.tsv"
#     log:
#         "logs/fastqc/{sample}.{read}.log"
#     threads: 1
#     resources:
#         mem="20G",
#         runtime="3h"
#     wrapper:
#         "v4.7.2/bio/fastqc"

rule fastqc:
    input:
        lookup(query="sample == '{sample}'", within=SAMPLES, cols="full_path"),
    output:
        html="qc/fastqc/{sample}.html",
        zip="qc/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    benchmark:
        "benchmarks/fastqc/{sample}.tsv"
    log:
        "logs/fastqc/{sample}.log"
    threads: 1
    resources:
        mem="10G",
        runtime="3h"
    wrapper:
        "v4.7.2/bio/fastqc"


rule samtools_stats:
    input:
        "mapped/{sample}.{type}.bam"
    params:
        extra=" ".join(["-t",config["reference"]["region_file"]]),
        region=""
    output:
        "qc/samtools-stats/{sample}.{type}.txt"
    resources:
        mem="10G",
        runtime="1h"
    benchmark:
        "benchmarks/samtools_stats/{sample}.{type}.tsv"
    log:
        "logs/samtools-stats/{sample}.{type}.log"
    wrapper:
        "v3.3.3/bio/samtools/stats"


rule picard_collect_hs_metrics:
    input:
        bam="mapped/{sample}.{type}.bam",
        reference="refs/genome.fasta",
        # Baits and targets should be given as interval lists. These can
        # be generated from bed files using picard BedToIntervalList.
        bait_intervals="refs/region.intervals",
        target_intervals="refs/region.intervals"
    output:
        "qc/hs_metrics/{sample}.{type}.txt"
    params:
        # Optional extra arguments. Here we reduce sample size
        # to reduce the runtime in our unit test.
        extra="--SAMPLE_SIZE 1000"
    benchmark:
        "benchmarks/picard_collect_hs_metrics/{sample}.{type}.tsv"
    log:
        "logs/picard/collect_hs_metrics/{sample}.{type}.log"
    resources:
        mem_mb=1024,
    wrapper:
        "v5.5.0/bio/picard/collecthsmetrics"
#        "v4.7.4-48-g99109e0/bio/picard/collecthsmetrics"

#rule multiqc_alignments:
#    shell:
#        """
#            multiqc {params} --force -o qc -n multiqc_alignments {input}
#        """

#rule multiqc_reads:
    #shell:
    #    """
    #       multiqc {params} -o qc -n multiqc_reads qc/fastqc/ >> {log} 2>&1
    #    """

rule multiqc:
    input:
        # expand("qc/fastqc/{sample}.{reads}_fastqc.zip", sample=SAMPLES.get("sample"), reads=["R1","R2"]),
        expand("qc/fastqc/{sample}_fastqc.zip", sample=SAMPLES.get("sample")),
        expand("qc/{ctype}/{sample}.{ftype}.txt", sample=SAMPLES.get("sample"), ctype=["samtools-stats","hs_metrics"], ftype=["woconsensus", "realigned"])
    output:
        "qc/multiqc.html",
        "qc/multiqc_data.zip",
    benchmark:
        "benchmarks/multiqc/multiqc.tsv"
    log:
        "logs/multiqc/multiqc.log"
    params:
        "--interactive --force --cl_config 'max_table_rows: 10000'"
    resources:
        mem="40G",
        runtime="2h"
    wrapper:
        "v5.5.2/bio/multiqc"
