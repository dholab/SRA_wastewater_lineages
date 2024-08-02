"""
TODO
"""

import os
import glob
import time
from dataclasses import dataclass, field
from loguru import logger
from .dereplicate import dereplicate
from typing import List


@dataclass(frozen=True)
class FetchBundle:
    """
    TODO
    """

    sra_acc: str
    bbmerge_path: str = "java -ea -Xmx8000m -Xms8000m bbmerge.sh"
    bbmap_path: str = "java -ea -Xmx8000m -Xms8000m bbmap.sh"
    ref: str = "SARS2.fasta"
    files: List[str] = field(init=False)

    def __post_init__(self):
        self.files = fetch_sra_accession(self)


def fetch_sra_accession(fetch_bundle: FetchBundle) -> List[str]:
    """
    TODO
    """
    logger.info(f"Beginning to fetch the accession {fetch_bundle.sra_acc} from SRA.")

    # run SRA prefetch
    os.system(f"prefetch {fetch_bundle.sra_acc}")
    os.system(f"fasterq-dump {fetch_bundle.sra_acc} --split-3")
    time.sleep(5)

    return glob.glob("*.fastq")


def handle_available_fastqs(available_fastqs: List[str], fetch_bundle: FetchBundle):
    """
    TODO
    """
    output_fasta = f"{fetch_bundle.sra_acc}.collapsed.fa"
    # make sure the file isn't orphaned
    if os.path.isfile(f"{fetch_bundle.sra_acc}_2.fastq"):
        logger.warning(
            f"""
            An orphaned R2 file was downloaded for {fetch_bundle.sra_acc}, {fetch_bundle.sra_acc}_2.fastq.
            The pipeline's current default behavior is to skip these cases.
            """
        )
        return

    # handle the accession based on whether it contains paired or unpaired reads
    if os.path.isfile(f"{fetch_bundle.sra_acc}_1.fastq") and os.path.isfile(
        f"{fetch_bundle.sra_acc}_2.fastq"
    ):
        logger.info(f"Processing {fetch_bundle.sra_acc} as paired reads.")

        # merge any overlapping reads, putting unmerged forward and reverse reads into
        # their own FASTQs
        os.system(
            f"{fetch_bundle.bbmerge_path} \
            qtrim=t \
            in1={fetch_bundle.sra_acc}_1.fastq in2={fetch_bundle.sra_acc}_2.fastq \
            out={fetch_bundle.sra_acc}.merge.fastq.gz \
            outu1={fetch_bundle.sra_acc}.unmerged_r1.fastq.gz \
            outu2={fetch_bundle.sra_acc}.unmerged_r2.fastq.gz"
        )

        # do away with the original FASTQs
        os.system(f"rm -f {fetch_bundle.sra_acc}_1.fastq")
        os.system(f"rm -f {fetch_bundle.sra_acc}_2.fastq")

        # Combine merged and unmerged reads into a single FASTQ
        os.system(
            f"cat {fetch_bundle.sra_acc}.merge.fastq.gz \
            {fetch_bundle.sra_acc}.unmerged_r1.fastq.gz \
            {fetch_bundle.sra_acc}.unmerged_r2.fastq.gz \
            > {fetch_bundle.sra_acc}.all.fastq.gz"
        )

        # clean up clean up everybody everywhere
        os.system(f"rm -f {fetch_bundle.sra_acc}.merge.fastq.gz")
        os.system(f"rm -f {fetch_bundle.sra_acc}.unmerged_r1.fastq.gz")
        os.system(f"rm -f {fetch_bundle.sra_acc}.unmerged_r2.fastq.gz")

        # if a mysterious third FASTQ came along, add that in too
        if os.path.isfile(fetch_bundle.sra_acc + ".fastq"):
            logger.info(
                "Combining all available reads from the three downloaded FASTQs."
            )
            os.system(
                f"cat {fetch_bundle.sra_acc}.fastq | gzip -c - >> {fetch_bundle.sra_acc}.all.fastq.gz"
            )
            os.system("rm -f " + fetch_bundle.sra_acc + ".fastq")

        # Dereplicate the reads
        logger.info("Dereplicating the reads")
        dereplicate(
            f"{fetch_bundle.sra_acc}.all.fastq.gz",
            f"{output_fasta}",
        )
        os.system(f"rm -f {fetch_bundle.sra_acc}.all.fastq.gz")

    elif os.path.isfile(f"{fetch_bundle.sra_acc}.fastq") or os.path.isfile(
        f"{fetch_bundle.sra_acc}_1.fastq"
    ):
        logger.info(f"Processing {fetch_bundle.sra_acc} as singleton reads.")
        dereplicate(f"{fetch_bundle.sra_acc}.fastq", f"{output_fasta}")
        os.system("rm -f " + fetch_bundle.sra_acc + ".fastq")

    elif os.path.isfile(f"{fetch_bundle.sra_acc}_1.fastq"):
        logger.info(f"Processing {fetch_bundle.sra_acc} as singleton reads.")
        dereplicate(f"{fetch_bundle.sra_acc}_1.fastq", f"{output_fasta}")
        os.system("rm -f " + fetch_bundle.sra_acc + ".fastq")
    else:
        logger.warning(
            f"Unsupported FASTQ name downloaded for {fetch_bundle.sra_acc}. Skipping"
        )
        return


def align_derep_reads(input_fasta: str, fetch_bundle: FetchBundle):
    """
    TODO
    """
    # run final mapping to the reference in preparation for haplotype calling
    logger.info("Mapping uncompressed file for {fetch_bundle.sra_acc}.")
    os.system(
        f"""
        minimap2 -a {fetch_bundle.ref} \
        {input_fasta} \
        --sam-hit-only --secondary=no \
        -o {fetch_bundle.sra_acc}.wg.sam
        """
    )
    os.system(f"rm -f {input_fasta}")
    logger.info(f"Finished fetching and deplicating {fetch_bundle.sra_acc}.")
