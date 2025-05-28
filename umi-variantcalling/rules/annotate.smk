rule annovar:
    input:
        "table/{sample}.csv"
    output:
        anno = temp("table/{sample}.anno.csv"),
        tmp = temp("table/{sample}.annotmp.csv")
    benchmark:
        "benchmarks/annovar/{sample}.tsv"
    log:
        "logs/annovar/{sample}.log"
    threads:
        8
    resources:
        runtime="30m",
        mem="20G"
    params:
        annovar = config['tools']['annovar'],
        anno = get_anno_params,
        sd = config["general"]['snakedir'],
        anno_info = f"scripts/shell/sed.cmd"
    shell:
        r"""
        tail -n +2 {input} > {output.tmp} && \
        {params.annovar}/table_annovar.pl {output.tmp} {params.anno} -thread {threads} && \
        sed --file={params.sd}/{params.anno_info} {output.tmp}.hg38_multianno.txt > {output.anno} && \
        rm {output.tmp}.hg38_multianno.txt
        """

# uncomment to add collums "EBScore" and "MultiAllelic", note that rule "add_ebfilter" and "umi-variantcalling/scripts/AddParameters.R" have to be adjusted too
# if config["general"]["control"]:
#     rule ebfilter:
#         input:
#             sample = lambda wildcards: input_bam[wildcards.sample],
#             vcf = "vardict/{sample}.vcf"
#         output:
#             vcf = temp("vardict/{sample}_EB.vcf"),
#             txt = temp("vardict/{sample}_EB.txt")
#         log:
#             "logs/EBFilter/{sample}.log"
#         threads:
#             4
#         resources:
#             time=get_time_3_1,
#             mem=get_mem_30_10
#         benchmark:
#             "benchmarks/ebfilter/{sample}.tsv"
#         conda:
#             "../env/EBFilter-env.yaml"
#         params:
#             normals = config['edit']['normals']
#         shell:
#             r"""
#             EBFilter -f vcf -t {threads} {input.vcf} {input.sample} {params.normals} {output.vcf}
#             bcftools query -f '[%EB]\n' {output.vcf} > {output.txt} 2>/dev/null
#             """
# else:
#     rule fake_ebfilter:
#         input:
#             vcf = "vardict/{sample}.vcf"
#         output:
#             txt = temp("vardict/{sample}_EB.txt")
#         log:
#             "logs/EBFilter/{sample}.log"
#         threads:
#             1
#         resources:
#             time=get_time_1_1
#         benchmark:
#             "benchmarks/fake_ebfilter/{sample}.tsv"
#         run:
#            with open(input.vcf, "r") as input_file:
#               vcf_lines = sum(1 for line in input_file if not line.startswith("#"))
           
#            with open(output.txt, "w") as output_file:
#               output_file.write("\n".join(["NaN"] * vcf_lines))

rule add_ebfilter:
    input:
        anno = "table/{sample}.anno.csv",
        #ebfilter = "vardict/{sample}_EB.txt" # uncomment if using EBFilter
    output:
        temp("table/{sample}.edit.csv")
    conda:
        "../env/Renv.yaml"
    benchmark:
        "benchmarks/add_ebfilter/{sample}.tsv"
    threads:
        1
    resources:
        runtime="30m"
    params:
        # wd = config["general"]["work_dir"] + "variantcalling/",
        sd = config["general"]['snakedir'],
        rs = f"scripts/AddParameters.R",
        candidate = config['edit']['candidate_list'],
        driver = config['edit']['driver_list'],
        CHIP = config['edit']['CHIP_list']
    shell:
        r"""
        Rscript {params.sd}/{params.rs} {input.anno} {output} {params.candidate} {params.driver} {params.CHIP} {input.ebfilter}
        """
        #Rscript {params.sd}/{params.rs} {params.wd}/{input.anno} {params.wd}/{output} {params.candidate} {params.driver} {params.CHIP} {params.wd}/{input.ebfilter}

rule primer3:
    input: "table/{sample}.edit.csv"
    output: temp("table/{sample}.edit.primer.csv")
    conda:
        "../env/primer3-env.yaml"
    threads: 1
    resources:
        runtime="30m",
        mem_mb=get_mem_3_1
    benchmark:
        "benchmarks/primer3/{sample}.tsv"
    params:
        genome_split = config["primer3"]["split"]
    script:
        "../scripts/primer3.py"


rule detect_HDR:
    input:
        filter_file = "table/{sample}.edit.csv",
        bam = "filterbam/{sample}.bam",
        index = "filterbam/{sample}.bai",
        pileup = "pileup/{sample}.pileup"
    output:
        HDR = temp("table/{sample}.edit.HDR.csv")
    benchmark:
        "benchmarks/detect_HDR/{sample}.tsv"
    conda:
        f"../env/HDR-env.yaml"
    threads: 1
    resources:
        time=get_time_3_1,
        mem=get_mem_40_10,
        mem_mb=get_mem_40_10
    params:
        min_sim = config['HDR']['min_similarity'],
        min_q = config['HDR']['min_q'],
        min_HDR_count = config['HDR']['min_HDR_count']
    script:
        "../scripts/HDR.py"


rule combine_annotations:
    input:
        anno_edit = "table/{sample}.edit.csv",
        primer = "table/{sample}.edit.primer.csv",
        hdr = "table/{sample}.edit.HDR.csv"
    output: 
        "filter/{sample}.edit.csv"
    benchmark:
        "benchmarks/combine_annotations/{sample}.tsv"
    conda:
        "../env/Renv.yaml"
    params: 
        wd = config["general"]["work_dir"] + "variantcalling/",
        sd = config["general"]['snakedir'],
        rs = f"scripts/CombineAnno.R"
    threads: 1
    resources:
        runtime="1h"
    shell:
        r"""
        Rscript {params.sd}/{params.rs} {input.anno_edit} {input.primer} {input.hdr} {output}
        """

# rule combine_samples:
#     input:
#         expand("filter/{sample}.edit.csv", sample = samples)
#     output: 
#         update("filter/variantcalls.csv")
#     conda:
#         "../env/Renv.yaml"
#     benchmark:
#         "benchmarks/combine_samples.tsv"
#     params: 
#         wd = config["general"]["work_dir"] + "variantcalling/",
#         sd = config["general"]['snakedir'],
#         rs = f"scripts/MergeCalls.R"
#     threads: 1
#     resources:
#         runtime="1h",
#         mem_mb=get_mem_30_10
#     shell:
#         r"""
#         Rscript {params.sd}/{params.rs} filter {output}
#         """
#         #Rscript {params.sd}/{params.rs} {params.wd}/filter {params.wd}/{output}

rule combine_samples:
    input:
        expand("filter/{sample}.edit.csv", sample = samples)
    output: 
        update("filter/variantcalls.csv")
    benchmark:
        "benchmarks/combine_samples.tsv"
    params: 
        wd = config["general"]["work_dir"] + "variantcalling/",
        sd = config["general"]['snakedir'],
        rs = f"scripts/MergeCalls.R"
    threads: 1
    resources:
        runtime="1h",
        mem="40G"
    shell:
        r"""
        FIRST_FILE=true

        # Process each .edit.csv file in the directory
        for FILE in {input}; do
            # Extract the sample name
            SAMPLE=$(basename "$FILE" .edit.csv)

            # Combine files
            if $FIRST_FILE; then
                # First file: copy header and content
                awk -v sample="$SAMPLE" 'BEGIN {{FS=OFS="\t"}} NR == 1 {{print "Sample", $0; next}} {{print sample, $0}}' "$FILE" > {output}
                FIRST_FILE=false
            else
                # Append only the content (skip header)
                awk -v sample="$SAMPLE" 'BEGIN {{FS=OFS="\t"}} NR == 1 {{print "Sample", $0; next}} {{print sample, $0}}' "$FILE" | tail -n +2 >> {output}
            fi
        done
        """



rule save:
    input:
        csv="filter/variantcalls.csv",
        vcf=expand("vardict/{sample}.vcf.gz", sample=samples),
    output:
        csv=config["longterm_storage"]["csv"],
        vcf=directory(config["longterm_storage"]["vcf"]),
        flag=touch("logs/save/varcall/done.txt")
    benchmark:
        "benchmarks/save/varcall.tsv"
    resources:
        runtime="12h"
    log:
        "logs/save/save_variantcalls.log"
    shell:
        """
        rsync -a {input.csv} {output.csv}
        rsync -a {input.vcf} {output.vcf} &> {log}
        """
