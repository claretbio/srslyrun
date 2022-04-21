rule fgbio_consensus:
    version:
        1
    message:
        "Consensus reads from UMIs with fgbio on {input}"
    input:
        "{base}.fgbio.bam"
    output:
        "{base}.fgbio.cons.bam"
    params:
        extra="-M 3 --tag=RX"
    log:
        "{base}-fgbio_consensus.log"
    wrapper:
        "v1.2.0/bio/fgbio/callmolecularconsensusreads"


rule fgbio_group:
    input:
        "{base}.rx.dupbam"
    output:
        bam="{base}.fgbio.bam",
        hist="{base}.fgbio.histo.tsv",
    params:
        extra="--s adjacency"
    log:
        "{base}-fgbio_group.log"
    wrapper:
        "v1.2.0/bio/fgbio/groupreadsbyumi"


rule index_rxdupbam:
    output:
        "{base}.rx.dupbam.bai",
    input:
        "{base}.rx.dupbam",
    wrapper:
        "v1.2.0/bio/samtools/index"


rule umi_tools_group:
    version:
        1
    message:
        "Correcting and grouping UMIs with UMI-Tools on {input}"
    output:
        "{base}.umi.dupbam",
        "{base}.umitools.tsv",
    input:
        "{base}.rx.dupbam",
        "{base}.rx.dupbam.bai",
    log:
        "{base}-umitools.log",
        "{base}-umitools.err"
    params:
        min_map_q=0,
        random_seed=42,
    shell:
        "umi_tools group"
        " --output-bam"
        " --stdin={input}"
        " --stdout={output[0]}"
        " --group-out={output[1]}"
        " --log={log[0]}"
        " --error={log[1]}"
        " --extract-umi-method=tag"
        " --umi-tag=RX"
        " --method=directional"
        " --paired"
        " --unmapped-reads=use"
        " --mapping-quality={params.min_map_q}"
        " --random-seed={params.random_seed}"


rule umi_srslyumi_bamtag:
    version:
        1
    message:
        "Moving UMI from fragment name to RX tag on {input}"
    output:
        temp("{base}.rx.dupbam"),
    input:
        "{base}.dupbam",
    log:
        "{base}.rxdupbam.log",
    shell:
        "srslyumi-bamtag --binary -o {output} --take-fragment 1 {input} 2> {log}"
