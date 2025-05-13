rule map_reads1a:
    input:
        unmapped=lookup(query="sample == '{sample}'", within=SAMPLES, cols="full_path"),
        #genome="refs/genome.fasta",
        #genomeindex="refs/genome.fasta.fai",
        #genomedict="refs/genome.dict",
        #bwaindex=multiext("refs/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        temp("mapped/{sample}.output.fastq"),
    conda:
        "../envs/picard.yaml"
    params:
        musage=config["picard"]["memoryusage"],
        #temp_file = temp("mapped/{sample}.temp.bam")
    group:
        "map_reads1"
    benchmark:
        "benchmarks/map_reads1a/{sample}.tsv"
    log:
        "logs/mapping/{sample}.fastq.log"
    threads:
        8
    resources:
        mem=get_mem_70_10,
        runtime="1h",
        mem_mib=7630
    shell:
        """
        picard {params.musage} SamToFastq -I {input.unmapped} -F {output} -INTERLEAVE true &> {log}
        """

rule map_reads1b: #bwa_mem2
    input:
        reads="mapped/{sample}.output.fastq",
        idx=multiext("refs/genome.fasta", ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123")
    output:
        temp("mapped/{sample}.bwaaligned.bam"),
    group:
        "map_reads1"
    benchmark:
        "benchmarks/map_reads1b_bwa2/{sample}.tsv"
    log:
        "logs/mapping/{sample}.bwaaligned_bwa2.log"
    threads:
        30
    resources:
        mem=get_mem_70_10,
        runtime="1h",
        mem_mib=7630
    params:
        extra=r"-p -R '@RG\tID:{sample}\tSM:{sample}'",
        sort="picard",  # Can be 'none', 'samtools', or 'picard'.
        sort_order="coordinate",  # Can be 'coordinate' (default) or 'queryname'.
        sort_extra="",  # Extra args for samtools/picard sorts.
    wrapper:
        "v5.5.0/bio/bwa-mem2/mem"

rule map_reads1c:
    input:
        unmapped=lookup(query="sample == '{sample}'", within=SAMPLES, cols="full_path"),
        genome="refs/genome.fasta",
        genomeindex="refs/genome.fasta.fai",
        genomedict="refs/genome.dict",
        bwaindex=multiext("refs/genome.fasta", ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123"),
        bwaaligned="mapped/{sample}.bwaaligned.bam"
    output:
        temp("mapped/{sample}.woconsensus1c.bam"),
    group:
        "map_reads1"
    conda:
        "../envs/mapping.yaml"
    params:
        musage=config["picard"]["memoryusage"],
        temp_file = temp("mapped/{sample}.temp.bam")
    benchmark:
        "benchmarks/map_reads1c/{sample}.tsv"
    log:
        "logs/mapping/{sample}.woconsensus1c.log"
    threads:
        8
    resources:
        mem=get_mem_70_10,
        runtime="1h",
        mem_mib=7630
    shell:
        """
        picard {params.musage} MergeBamAlignment \
        UNMAPPED={input.unmapped} ALIGNED={input.bwaaligned} O={output} R={input.genome} \
        SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT &> {log}
        """

rule replace_rg1a:
    input:
        "mapped/{sample}.woconsensus1c.bam"
    output:
        temp("mapped/{sample}.woconsensus.woRG.bam")
    conda:
        "../envs/samtools.yaml"
    group:
        "map_reads1"
    resources:
        runtime="30m",
#        mem="2G"
    log:
        "logs/replace_rg1a/{sample}.log"
    benchmark:
        "benchmarks/replace_rg1a/{sample}.tsv"
    shell:
        """
        samtools view -h {input} | grep -v "^@RG" | sed "s/\tRG:Z:[^\t]*//" | samtools view -bo {output} &> {log}
        """

rule replace_rg1b:
    input:
        "mapped/{sample}.woconsensus.woRG.bam",
    output:
        temp("mapped/{sample}.woconsensus.bam"),
    group:
        "map_reads1"
    log:
        "logs/replace_rg1b/{sample}.log",
    benchmark:
        "benchmarks/replace_rg1b/{sample}.tsv"
    params:
        extra="--RGLB illumina --RGPL illumina --RGPU {sample} --RGSM {sample} --RGID {sample}",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
 #       mem="2G",
        runtime="30m"
    wrapper:
        "v5.5.0/bio/picard/addorreplacereadgroups"

rule GroupReads:
    input:
        "mapped/{sample}.woconsensus.bam"
    output:
        bam=temp("unmapped/{sample}.groupedumi.bam"),
        hist="metrices/fgbio/{sample}.groupedumi.histo.tsv",
    params:
        extra=config["fgbio"]["groupreads"]
    benchmark:
        "benchmarks/GroupReads/{sample}.tsv"
    resources:
        mem="30G",
        runtime="2h"
    log:
        "logs/fgbio/group_reads/{sample}.log"
    wrapper:
        "v2.4.0/bio/fgbio/groupreadsbyumi"

rule ConsensusReads:
    input:
        "unmapped/{sample}.groupedumi.bam"
    output:
        temp("unmapped/{sample}.consensusreads.bam")
        #temporary("unmapped/{sample}.consensusreads.bam")
    params:
        extra=config["fgbio"]["callconsensus"]
    benchmark:
        "benchmarks/ConsensusReads/{sample}.tsv"
    resources:
        mem="10G",
        runtime="2h"
    log:
        "logs/fgbio/consensus_reads/{sample}.log"
    wrapper:
        "v2.4.0/bio/fgbio/callmolecularconsensusreads"

####
### Common errors:
# "Error in writing fastq file /dev/stdout"
# means the pipeline collapsed. Probably an error in MergeBamAlignment

#rule map_reads2:
#    input:
#        unmapped="unmapped/{sample}.consensusreads_qsorted.bam",
#        genome="refs/genome.fasta",
#        genomeindex="refs/genome.fasta.fai",
#        genomedict="refs/genome.dict",
#        bwaindex=multiext("refs/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
#    output:
#        temp("mapped/{sample}.consensusreads.bam")
#    conda:
#        "../envs/mapping.yaml"
#    params:
#        musage=config["picard"]["memoryusage"],
#    benchmark:
#        "benchmarks/map_reads2/{sample}.tsv"
#    log:
#        "logs/mapping/{sample}.consensusreads.log"
#    threads:
#        8
#    resources:
#        mem=get_mem_100_50,
#        runtime="5h"
#    shell:
#        r"""
#        picard {params.musage} SamToFastq -I {input.unmapped} -F /dev/stdout -INTERLEAVE true \
#        | bwa mem -p -t 8 {input.genome} /dev/stdin \
#        | picard {params.musage} MergeBamAlignment \
#        UNMAPPED={input.unmapped} ALIGNED=/dev/stdin O={output} R={input.genome} \
#        SO=queryname CREATE_INDEX=true ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT &> {log}
#        """

rule map_reads2a:
    input:
        unmapped="unmapped/{sample}.consensusreads_qsorted.bam",
    output:
        temp("mapped/{sample}.output2.fastq"),
    conda:
        "../envs/picard.yaml"
    params:
        musage=config["picard"]["memoryusage"],
    group:
        "map_reads2"
    benchmark:
        "benchmarks/map_reads2a/{sample}.tsv"
    log:
        "logs/mapping/{sample}.fastq2.log"
    threads:
        8
    resources:
        mem=get_mem_70_10,
        runtime="1h",
        mem_mib=7630
    shell:
        """
        picard {params.musage} SamToFastq -I {input.unmapped} -F {output} -INTERLEAVE true &> {log}
        """

rule map_reads2b: #bwa_mem2
    input:
        reads="mapped/{sample}.output2.fastq",
        idx=multiext("refs/genome.fasta", ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123")
    output:
        temp("mapped/{sample}.bwaaligned2.bam"),
    group:
        "map_reads2"
    benchmark:
        "benchmarks/map_reads2b_bwa2/{sample}.tsv"
    log:
        "logs/mapping/{sample}.bwaaligned2_bwa2.log"
    threads:
        30
    resources:
        mem=get_mem_70_10,
        runtime="1h",
        mem_mib=7630
    params:
        extra=r"-p -R '@RG\tID:{sample}\tSM:{sample}\tLB:targeted\tPU:flowcell\tPL:ILLUMINA'",
        sort="picard",  # Can be 'none', 'samtools', or 'picard'.
        sort_order="queryname",  # Can be 'coordinate' (default) or 'queryname'.
        sort_extra="",  # Extra args for samtools/picard sorts.
    wrapper:
        "v5.5.0/bio/bwa-mem2/mem"

rule map_reads2c:
    input:
        unmapped="unmapped/{sample}.consensusreads_qsorted.bam",
        genome="refs/genome.fasta",
        genomeindex="refs/genome.fasta.fai",
        genomedict="refs/genome.dict",
        bwaindex=multiext("refs/genome.fasta", ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123"),
        bwaaligned="mapped/{sample}.bwaaligned2.bam"
    output:
        temp("mapped/{sample}.consensusreads.bam"),
    conda:
        "../envs/mapping.yaml"
    params:
        musage=config["picard"]["memoryusage"],
        temp_file = temp("mapped/{sample}.temp.bam")
    group:
        "map_reads2"
    benchmark:
        "benchmarks/map_reads2c/{sample}.tsv"
    log:
        "logs/mapping/{sample}.consensusreads.log"
    threads:
        8
    resources:
        mem=get_mem_70_10,
        runtime="1h",
        mem_mib=7630
    shell:
        """
        picard {params.musage} MergeBamAlignment \
        UNMAPPED={input.unmapped} ALIGNED={input.bwaaligned} O={output} R={input.genome} \
        SO=queryname CREATE_INDEX=true ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT &> {log}
        """


rule FilterConsensusReads:
    input:
        "mapped/{sample}.consensusreads.bam"
    output:
        temp("mapped/{sample}.filtered.bam")
    conda: "../envs/fgbio.yaml"
    params:
        extra=config["fgbio"]["fextra"],
        min_base_quality=config["fgbio"]["fminq"],
        min_reads=[3],
        ref="refs/genome.fasta",
        musage=config["picard"]["memoryusage"],
    benchmark:
        "benchmarks/FilterConsensusReads/{sample}.tsv"
    log:
        "logs/fgbio/filterconsensusreads/{sample}.log"
    threads: 1
    resources:
        mem="50G",
        runtime="1h"
    shell:
        r"""
        fgbio {params.musage} FilterConsensusReads \
         -i {input} \
         -o {output} \
         -r {params.ref} \
         --min-reads {params.min_reads} \
         --min-base-quality {params.min_base_quality} \
         &> {log}
         """

rule replace_rg2a:
    input:
        "mapped/{sample}.filtered.bam"
    output:
        temp("mapped/{sample}.filtered.woRG.bam")
    conda:
        "../envs/samtools.yaml"
    resources:
        runtime="30m",
        mem="2G"
    log:
        "logs/replace_rg2a/{sample}.log"
    benchmark:
        "benchmarks/replace_rg2a/{sample}.tsv"
    shell:
        """
        samtools view -h {input} | grep -v "^@RG" | sed "s/\tRG:Z:[^\t]*//" | samtools view -bo {output} &> {log}
        """

rule replace_rg2b:
    input:
        "mapped/{sample}.filtered.woRG.bam",
    output:
        temp("fixed-rg/{sample}.filtered.bam"),
    log:
        "logs/replace_rg2b/{sample}.log",
    benchmark:
        "benchmarks/replace_rg2b/{sample}.tsv"
    params:
        extra="--RGLB illumina --RGPL illumina --RGPU {sample} --RGSM {sample} --RGID {sample}",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem="2G",
        runtime="30m"
    wrapper:
        "v5.5.0/bio/picard/addorreplacereadgroups"


rule realignertargetcreator:
    input:
        bam="mapped/{sample}.filtered_csorted.bam",
        bam_index="mapped/{sample}.filtered_csorted.bam.bai",
        bed=config["reference"]["region_file"],
        ref="refs/genome.fasta",
        known="refs/known_indels.vcf.gz"
    output:
        temp("realigned/{sample}.intervals")
    conda: "../envs/gatk3.yaml"
    benchmark:
        "benchmarks/realignertargetcreator/{sample}.tsv"
    log:
        "logs/gatk/realignertargetcreator/{sample}.log"
    params:
        extra="",  # optional
        java_opts=config["picard"]["memoryusage"],
    threads: 10
    resources:
        mem="70G",
        runtime="2h"
    shell:
        r"""
        gatk3 {params.java_opts} -T RealignerTargetCreator {params.extra} -nt {threads} -I {input.bam} -R {input.ref} -known {input.known} -L {input.bed} -o {output} &> {log}
        """


rule indelrealigner:
    input:
        bam="mapped/{sample}.filtered_csorted.bam",
        bam_index="mapped/{sample}.filtered_csorted.bam.bai",
        bed=config["reference"]["region_file"],
        ref="refs/genome.fasta",
        known="refs/known_indels.vcf.gz",
        target_intervals="realigned/{sample}.intervals"
    output:
        "mapped/{sample}.realigned.bam"
    conda: "../envs/gatk3.yaml"
    benchmark:
        "benchmarks/indelrealigner/{sample}.tsv"
    log:
        "logs/gatk3/indelrealigner/{sample}.log"
    params:
        extra="",  # optional
        java_opts=config["picard"]["memoryusage"], # optional
    threads: 5
    resources:
        mem="70G",
        runtime="5h"
    shell:
        """
        gatk3 {params.java_opts} -T IndelRealigner {params.extra} -I {input.bam} -R {input.ref} -known {input.known} -L {input.bed} --targetIntervals {input.target_intervals} -o {output}
        """
