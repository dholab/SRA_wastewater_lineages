#!/bin/env python3

import os
import time
from argparse import ArgumentParser


def parse_process_arrays_args(parser: ArgumentParser):
    """Parses the python script arguments from bash and makes sure files/inputs are valid"""
    parser.add_argument(
        "--SRA",
        type=str,
        help="SRA number must be case-sensitive perfect match",
        required=True,
    )


def get_process_arrays_args():
    """Inputs arguments from bash
    Gets the arguments, checks requirements, returns a dictionary of arguments
    Return: args - Arguments as a dictionary
    """
    parser = ArgumentParser()
    parse_process_arrays_args(parser)
    args = parser.parse_args()
    return args


args = get_process_arrays_args()
SRA = args.SRA
BBMERGE_PATH = "java -ea -Xmx8000m -Xms8000m bbmerge.sh"
BBMAP_PATH = "java -ea -Xmx8000m -Xms8000m bbmap.sh"
DEREP_PATH = "./derep.py"
SARS2_REF = "SARS2.fasta"


def fetch(SRA_ID):
    # removed if already exists be cause it can't exist in the chtc image?
    print(SRA_ID)
    print(time.ctime(time.time()))
    os.system("prefetch " + SRA_ID)
    os.system("fasterq-dump " + SRA_ID + " --split-3")
    time.sleep(5)

    if os.path.isfile(SRA_ID + "_1.fastq") and os.path.isfile(SRA_ID + "_2.fastq"):
        print("--paired reads--")
        os.system(
            f"{BBMERGE_PATH} qtrim=t in1={SRA_ID}_1.fastq in2={SRA_ID}_2.fastq  out={SRA_ID}.merge.fq outu1={SRA_ID}.un1.fq outu2={SRA_ID}.un2.fq"
        )
        os.system("rm -f " + SRA_ID + "_1.fastq")
        os.system("rm -f " + SRA_ID + "_2.fastq")
        print("combining merged with unique")
        os.system(
            f"cat {SRA_ID}.merge.fq {SRA_ID}.un1.fq {SRA_ID}.un2.fq > {SRA_ID}.all.fq"
        )
        os.system("rm -f " + SRA_ID + ".merge.fq")
        os.system("rm -f " + SRA_ID + ".un1.fq")
        os.system("rm -f " + SRA_ID + ".un2.fq")
        if os.path.isfile(SRA_ID + ".fastq"):
            print("combining merged with unique fastq to all fastq")
            os.system(f"cat {SRA_ID}.fastq >> {SRA_ID}.all.fq")
            os.system("rm -f " + SRA_ID + ".fastq")
        print("Dereplicating the reads")
        os.system(f"python {DEREP_PATH} {SRA_ID}.all.fq {SRA_ID}.collapsed.fa 1")
        os.system("rm -f " + SRA_ID + ".all.fq")
    elif os.path.isfile(SRA_ID + ".fastq"):
        print("--singleton reads--")
        print("Dereplicating the reads")
        os.system(f"python {DEREP_PATH} {SRA_ID}.fastq {SRA_ID}.collapsed.fa 1")
        os.system("rm -f " + SRA_ID + ".fastq")
    elif os.path.isfile(SRA_ID + "_1.fastq"):
        print("singleton reads")
        print("Dereplicating the reads")
        os.system(f"python {DEREP_PATH} {SRA_ID}_1.fastq {SRA_ID}.collapsed.fa 1")
        os.system("rm -f " + SRA_ID + "_1.fastq")
    elif os.path.isfile(SRA_ID + "_2.fastq"):
        print("------------------------------------------ ")
        print("------------------------------------------ ")
        print("------------------------------------------ ")
        print("Single pair orphanned")
        print("Single ")
        print("Single ")
        print("------------------------------------------ ")
        print("------------------------------------------ ")
        print("------------------------------------------ ")

    if os.path.isfile(SRA_ID + ".collapsed.fa"):
        print("mapping uncompressed file")
        os.system(
            "minimap2 -a "
            + SARS2_REF
            + " "
            + SRA_ID
            + ".collapsed.fa --sam-hit-only --secondary=no -o "
            + SRA_ID
            + ".SARS2.wg.sam"
        )
        os.system("rm -f " + SRA_ID + ".collapsed.fa")

    elif os.path.isfile(SRA_ID + ".collapsed.fa.gz"):
        print("mapping compressed file")
        os.system(
            "minimap2 -a "
            + SARS2_REF
            + " "
            + SRA_ID
            + ".collapsed.fa.gz --sam-hit-only --secondary=no -o "
            + SRA_ID
            + ".SARS2.wg.sam"
        )
        os.system("rm -f " + SRA_ID + ".collapsed.fa.gz")
    print(SRA_ID + " done")


fetch(SRA)
