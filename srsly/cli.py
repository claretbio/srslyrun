#!/usr/bin/env python3
"""
Read in a list of libraries and run basic informatic analysis
"""


import argparse
import glob
import json
import multiprocessing
import re
import os
import sys
import gzip
import yaml
import snakemake
import tarfile
import importlib.resources


def read_list(fn):
    with open(fn) as f:
        lines = [x.rstrip("\r\n") for x in f]
    return lines


def read_csv(fn, field_delim=",", line_end="\r\n"):
    lines = []
    with open(fn) as f:
        for line in f:
            lines.append(line.rstrip(line_end).split(field_delim))
    return lines


class NoReferenceSpecified(Exception):
    pass


class NoLibrariesSpecified(Exception):
    pass


class BothLibfileAndCLILibraries(Exception):
    pass


class SampleSheetMissingError(Exception):
    def __init__(self, filename):
        self.filename = filename


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


def get_sampledict_from_cli_args(args):
    # TODO: gather additional columns with per-sample configuration, such as
    # enrichment panel, UMI, additional analysis modules, etc.
    if args.libraries:
        samples = [re.sub("-[0-9]{6}$", "", x) for x in args.libraries]
        return dict(zip(args.libraries, samples))
    if args.libfile:
        libraries = open(args.libfile).read().splitlines()
        samples = [re.sub("-[0-9]{6}$", "", x) for x in libraries]
        return dict(zip(libraries, samples))
    samplesheet_fn = args.samplesheet
    fastqdir = args.fastqdir or get_absolute_directory(args.samplesheet)
    seqdate = args.seqdate or parse_seq_date(fastqdir)
    try:
        return read_samplesheet(samplesheet_fn, seqdate)
    except FileNotFoundError:
        raise SampleSheetMissingError(samplesheet_fn)


def result_dir_for(results_dir, samplename):
    ## TODO: if breaking out samples by pipline (SR/APN), insert logic here
    os.path.join(results_dir, samplename)


###
### Run Pipeline on Samples
###
def get_samples(args):
    if args.libraries:
        samples = [re.sub("-[0-9]{6}$", "", x) for x in args.libraries]
        return dict(zip(args.libraries, samples))
    if args.libfile:
        libraries = open(args.libfile).read().splitlines()
        samples = [re.sub("-[0-9]{6}$", "", x) for x in libraries]
        return dict(zip(libraries, samples))


def run_sample_pipeline(args):
    if args.libfile and args.libraries:
        raise BothLibfileAndCLILibraries()
    libs = args.libfile or args.libraries
    if not libs:
        raise NoLibrariesSpecified()

    if not args.reference:
        raise NoReferenceSpecified()

    if not os.path.isdir('{}/workflow'.format(args.resultsdir)):
        with importlib.resources.path('srsly', 'workflow.tar.gz') as tarball:
            tar=tarfile.open(tarball)
            tar.extractall(path=args.resultsdir)


    configfilename = write_run_config(args)
    libraries = get_samples(args)
    customconfig = {"libraries": ",".join(libraries.keys())}
    if args.resultsdir:
        customconfig["resultsdir"] = args.resultsdir
    if args.umi:
        customconfig["umi"] = "true"
    else:
        customconfig["umi"] = "false"

    snakemake.snakemake(
        snakefile=sample_pipeline_snakefile(args),
        workdir=args.resultsdir,
        cores=args.cores or multiprocessing.cpu_count(),
        dryrun=args.dryrun,
        printshellcmds=True,
        use_conda=True,
        printreason=True,
        configfiles=[configfilename],
        config=customconfig,
        cleanup_metadata=args.cleanup,
        resources={"mem_mb": 150000},
        printdag=args.dag,
        debug_dag=args.debug_dag,
        unlock=args.unlock,
        force_incomplete=args.rerun_incomplete
    )


    

def repo_dir(args):
    return os.path.realpath(args.resultsdir)


def repo_config_filename(args):
    return os.path.join(repo_dir(args), "workflow", "config", "config.json")


def sample_pipeline_snakefile(args):
    return os.path.join(repo_dir(args), "workflow", "Snakefile")


def snakemake_config_path(args):
    return os.path.join(repo_dir(args), "workflow", "config", "config.json")


def get_absolute_directory():
    return os.path.dirname(os.path.abspath())


def commalist(x, sep=","):
    return x.split(",")


def listfile(fn):
    f = open(fn)
    l = [x.strip() for x in f]
    f.close()
    return l


def results_dir(resultsdir, config):
    if resultsdir:
        return resultsdir
    with open(config) as f:
        cfg = json.load(f)
    return cfg["resultsdir"]

def write_run_config(args):
    configdict = {}
    configdict['reference'] = os.path.realpath(args.reference)
    configdict['indir'] = os.path.realpath(args.fastqdir)
    configdict['resultsdir'] = os.path.realpath(args.resultsdir)
    configdict['umi'] = 'true' if args.umi else 'false'
    configfilename = '{}/workflow/config/config.json'.format(configdict['resultsdir'])
    with open(configfilename, 'w') as f:
        json.dump(configdict, f) 
    return configfilename

def add_sample_list_args(subparser):
    subparser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing config files"
    )
    subparser.add_argument(
        "--fastqdir",
        default=os.getcwd(),
        help="directory with raw FASTQ files. Defaults to the current "
        "working directory ",
    )
    subparser.add_argument(
        "--reference",
        help="reference genome for samples",
    )
    subparser.add_argument(
        "--libraries",
        type=commalist,
        help="Instead SampleSheet.csv, use this list of comma "
        "separated libname-seqdate libraries",
    )
    subparser.add_argument(
        "--libfile",
        help="Instead of listing on the command line, use this file "
        "containing one library identifier per line",
    )
    subparser.add_argument(
        "--resultsdir",
        default=os.getcwd(),
        help="specify output directory",
    )
    subparser.add_argument(
        "--umi",
        action="store_true",
        default=False,
        help="run the UMI aware version of analysis"
    )
    return subparser


def add_snakemake_args(subparser):
    subparser.add_argument(
        "--cleanup",
        type=read_list,
        help="Recover from an interrupted run on these files (e.g. power outage, Ctrl-C). The filelist should have been output from a prior invocation of runsamples, you'll have to copy and paste.",
    )

    subparser.add_argument(
        "--dryrun",
        action="store_true",
        help="Don't run any commands, just print what will be run",
    )
    subparser.add_argument(
        "--cores",
        type=int,
        help="Number of compute cores to use (default: all)",
    )
    subparser.add_argument(
        "--dag",
        default=False,
        action="store_true",
        help="create a dag from the rules to be run (Default=False)",
    )
    subparser.add_argument(
        "--debug_dag",
        default=False,
        action="store_true",
        help="print candidate and selected jobs and their wildcards to aid in debugging (Default=False)",
    )
    subparser.add_argument(
        "--unlock",
        default=False,
        action="store_true",
        help="unlock locked working directory after run interruption",
    )
    subparser.add_argument(
            "--rerun-incomplete",
            default=False,
            action="store_true",
            help="rerun jobs left incomplete by a run interruption",
            )
    return subparser


def argparser():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.set_defaults(run_command=None)
    sp = ap.add_subparsers(title="Operations", dest="operation")

    ## Run samples sub-command
    rs = sp.add_parser("runsamples", help="run the samples in a flowcell directory")
    rs.set_defaults(run_command=run_sample_pipeline)
    rs = add_sample_list_args(rs)
    rs = add_snakemake_args(rs)

    return ap


def exe_in_path(program, envvar="PATH"):
    import os

    for path in os.environ["PATH"].split(os.pathsep):
        exe_file = os.path.join(path, program)
        if os.path.isfile(exe_file) and os.access(exe_file, os.X_OK):
            return exe_file
    return None


def main():
    ap = argparser()
    a = ap.parse_args()

    if not exe_in_path("conda"):
        print(
            "Could not find `conda' in the PATH environment variable. Please install conda to continue"
        )
        return 7

    if a.run_command is None:
        ap.print_help()
        return 0
    try:
        run_command = a.run_command
        del a.operation
        del a.run_command
        run_command(a)
    except FastqMissingError as e:
        msg = "ERROR: couldn't find a FASTQ file in '{}' for sample: {}"
        print(msg.format(e.indir, e.sample))
        return 2
    except TooManyMatchingFastqError as e:
        msg = "ERROR: Found too many FASTQ for {} in '{}'"
        print(msg.format(e.sample, e.indir))
        for f in e.files:
            print("\t" + f)
        return 3
    except NoLibrariesSpecified:
        print("Neither --libraries nor --libfile were specified, one is required")
        return 4
    except BothLibfileAndCLILibraries:
        print("ERROR: both --libraries and --libfile were used, only one is allowed")
        return 5
    except NoReferenceSpecified:
        print("ERROR: no reference genome specified, this is required")
        return 6

if __name__ == "__main__":
    sys.exit(main())
