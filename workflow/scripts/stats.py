#!/usr/bin/env python3
import argparse
import sys


def read_cutadapt(fn):
    r = {}
    with open(fn) as f:
        for line in f:
            if line.startswith("Total read pairs"):
                key, _, val = line.rstrip().partition(":")
                r[key] = val.strip()
            if line.startswith("Pairs"):
                key, _, val = line.rstrip().partition(":")
                r[key] = val.strip()
    return r


def read_picard_metrics_file(fn):
    conv = {"LIBRARY": str, "PERCENT_DUPLICATION": float}
    with open(fn) as f:
        for line in f:
            if line.startswith("## METRICS CLASS"):
                break
        header = f.readline().rstrip().split("\t")
        values = f.readline().rstrip().split("\t")
        return dict([(h, conv.get(h, int)(v)) for h, v in zip(header, values)])


def read_flagstats(fn):
    r = {}
    with open(fn) as f:
        for line in f:
            count, _, _, *description = line.rstrip().split(" ")
            key = " ".join(description)
            if key.startswith("properly paired"):
                key = "properly paired"
            elif key.startswith("mapped"):
                key = "mapped"
            elif key.startswith("singletons"):
                key = "singletons"
            r[key] = int(count)
    r["total"] = (
        r["in total (QC-passed reads + QC-failed reads)"]
        - r["secondary"]
        - r["supplementary"]
    )
    r["fraction_pairs_mapped"] = (
        r["with itself and mate mapped"] + 2 * r["singletons"]
    ) / float(r["total"])
    r["fraction_chimeras"] = r["with mate mapped to a different chr"] / float(
        r["total"]
    )
    return r


def main(cutadaptfn, flagstatfn, picard_dupsfile, out):
    cutadapt = read_cutadapt(cutadaptfn)
    flagstat = read_flagstats(flagstatfn)
    picard = read_picard_metrics_file(picard_dupsfile)

    o = [
        ("pairs", cutadapt["Total read pairs processed"]),
        ("pairs_too_short", cutadapt["Pairs that were too short"]),
        ("pairs_kept", cutadapt["Pairs written (passing filters)"]),
        ("pairs_mapped", flagstat["with itself and mate mapped"] / 2),
        ("singletons_mapped", flagstat["singletons"]),
        (
            "fragments_mapped",
            flagstat["singletons"] + flagstat["with itself and mate mapped"] // 2,
        ),
        ("percent_bam_mapped", flagstat["fraction_pairs_mapped"] * 100),
        (
            "percent_fastq_mapped",
            (
                flagstat["total"]
                * flagstat["fraction_pairs_mapped"]
                / 2
                / int(cutadapt["Total read pairs processed"].replace(',', ''))
            )
            * 100,
        ),
        ("proper_pairs", flagstat["properly paired"] / 2),
        (
            "percent_bam_proper_pairs",
            (flagstat["properly paired"] / float(flagstat["total"])) * 100,
        ),
        (
            "picard_dups",
            picard["READ_PAIR_DUPLICATES"] + picard["UNPAIRED_READ_DUPLICATES"],
        ),
        ("picard_percent_dup", (picard["PERCENT_DUPLICATION"]) * 100),
        ("picard_est_libsize", picard.get("ESTIMATED_LIBRARY_SIZE", -1)),
    ]
    for key, val in o:
        out.write("{}\t{}\n".format(key, val))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "-o",
        default=sys.stdout,
        help="output file (stdout)",
        type=argparse.FileType("w"),
    )
    ap.add_argument("flagstat", help="samtools flagstat output")
    ap.add_argument("picard_dupmetrics", help="Picard tools duplicate metrics")
    ap.add_argument("cutadapt_paired", help="cutadapt qc fil")

    a = ap.parse_args()
    main(a.cutadapt_paired, a.flagstat, a.picard_dupmetrics,a.o)

