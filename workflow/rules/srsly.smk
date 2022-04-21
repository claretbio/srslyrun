import glob


rule stats:
    version:
       "1"
    message:
       "Collating stats to {output}"
    input:
        metrics="{base}.picardmetrics.txt",
        flagstat="{base}.bam.flagstat",
        cutadapt="{base}.qc.txt"
    output:
        "{base}.stats"
    shell:
        "workflow/scripts/stats.py -o {output}"
        " {input.flagstat}"
        " {input.metrics}"
        " {input.cutadapt}"


rule flagstat:
    version:
       "1"
    message:
       "Running samtools flagstat on {input}"
    input:
        "{base}.bam",
    output:
        "{base}.bam.flagstat",
    log: 
        "{base}-flagstat.log"
    wrapper:
        "v1.3.1/bio/samtools/flagstat"


rule plot_mol_len:
    version:
         "1"
    output:
        "{base}.insert_length_dist.pdf"
    input:
        "{base}.mol_lens.txt"
    #params:
    #    sample=lambda wildcards: wildcards.sample
    shell:
        "Rscript workflow/scripts/plot_mol_len.R {input} $(basename {wildcards.base}) {output}" 


rule write_mol_len:
    version:
        "1"
    output:
        temp("{base}.mol_lens.txt")
    input:
        bam="{base}.mml.bam",
    shell:
        "python3 workflow/scripts/write_mol_len.py {input.bam} {output}"


rule filter_mol_len:
    version:
        "1"
    output:
        bam=temp("{base}.mml.bam"),
        #bai=temp("{base}.mml.bam.bai")
    input:
        "{base}.bam"
    params:
        extra="-f64 -F3840 -q20",
        samtools_opts="--write-index"
    wrapper:
        "v1.3.1/bio/samtools/view"


rule bam_index:
    version:
        "1"
    message:
        "Indexing bam {input}"
    output:
        "{base}.bam.bai",
    input:
        "{base}.bam",
    log:
        "{base}-index.log",
    wrapper:
        "v1.3.1/bio/samtools/index"


def get_params(umi):
    if umi == 'true':
        return "--BARCODE_TAG BX --REMOVE_DUPLICATES true"
    else:
        return "--REMOVE_DUPLICATES true"


rule mark_duplicates:
    version:
        "1"
    message:
        "Marking duplicates with Picard for {input}"
    input:
        "{base}.dupbam",
    output:
        bam="{base}.bam",
        metrics="{base}.picardmetrics.txt",
    log:
        "{base}-picard.log",
    params:
        extra=get_params(config['umi'])
    resources:
        mem_mb=1024,
    wrapper:
        "v1.3.1/bio/picard/markduplicates"


if config['reference'].endswith('.fa'):
    ref = config['reference'].replace('.fa', '')
elif config['reference'].endswith('.fasta'):
    ref = config['reference'].replace('.fasta', '')


rule bwa_mem:
    version:
        "1"
    message:
        "Aligning with BWA MEM and sorting with samtools"
    output:
        temp("{base}.dupbam"),
    input:
        reads=['{base}_R1.fastq',
               '{base}_R2.fastq'],
        #idx="{reference}".format(**config)
        idx=multiext(ref, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    params:
        sorting="samtools",
    threads: 8
    log:
        "{base}-mem.log",
    wrapper:
        "v1.3.1/bio/bwa/mem"

class FastqMissingError(Exception):
    def __init__(self, read, base, rawdir):
        self.read = read
        self.base = base
        self.rawdir = rawdir


class TooManyMatchingFastqError(Exception):
    def __init__(self, read, base, rawdir, files):
        self.read = read
        self.base = base
        self.rawdir = rawdir
        self.files = files


def find_a_fastq(indir, sample, read):
    pattern = "{}/{}*{}*.fastq.gz".format(indir, sample, read)
    reads = glob.glob(pattern)

    if len(reads) == 0:
        raise FastqMissingError(indir, sample, read)
    elif len(reads) > 1:
        raise TooManyMatchingFastqError(indir, sample, read, reads)

    return reads[0]


def find_fastqs(sample):
    r1 = find_a_fastq(config["indir"], sample, "R1")
    r2 = find_a_fastq(config["indir"], sample, "R2")
    return [r1, r2]


rule cutadapt:
    version:
        "1"
    message:
        "Trimming adapters with cutadapt"
    output:
        fastq1="{resultsdir}/{{sample}}/{{sample}}_R1.fastq".format(**config),
        fastq2="{resultsdir}/{{sample}}/{{sample}}_R2.fastq".format(**config),
        qc="{resultsdir}/{{sample}}/{{sample}}.qc.txt".format(**config)
    input:
        lambda wildcards: find_fastqs(wildcards.sample)
    params:
        adapters="-a AGATCGGAAGAGCACACGTCTGAA -g AGATCGGAAGAGCGTCGTGTAGGG -A AGATCGGAAGAGCACACGTCTGAA -G AGATCGGAAGAGCGTCGTGTAGGG",
        extra="--minimum-length 30"
    log:
       "{resultsdir}/{{sample}}/{{sample}}-cutadapt.log".format(**config)
    wrapper:
        "v1.3.1/bio/cutadapt/pe"


#if config['reference'].endswith('.fa'):
#    ref = config['reference'].replace('.fa', '')
#elif config['reference'].endswith('.fasta'):
#    ref = config['reference'].replace('.fasta', '')


rule bwa_index:
    input:
        "{reference}".format(**config),
    output:
        idx=multiext(ref, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa-index.log",
    params:
        algorithm="bwtsw",
    wrapper:
        "v1.3.1/bio/bwa/index"
