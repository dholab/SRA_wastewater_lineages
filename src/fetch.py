"""
TODO
"""

from __future__ import annotations

import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Self

from loguru import logger

from .dereplicate import dereplicate


@dataclass
class FetchBundle:
    """
    TODO
    """

    sra_acc: str
    bbmerge_path: str = "java -ea -Xmx8000m -Xms8000m bbmerge.sh"
    bbmap_path: str = "java -ea -Xmx8000m -Xms8000m bbmap.sh"
    ref: str = "SARS2.fasta"
    files: list[Path] = field(init=False)

    def __post_init__(self) -> Self:
        self.files = fetch_sra_accession(self)
        return self


def fetch_sra_accession(fetch_bundle: FetchBundle) -> list[Path]:
    """
    TODO
    """
    logger.info(
        f"Beginning to fetch the accession {fetch_bundle.sra_acc} from SRA.",
    )

    # run SRA prefetch
    subprocess.run(["prefetch", fetch_bundle.sra_acc], check=True)
    subprocess.run(["fasterq-dump", fetch_bundle.sra_acc, "--split-3"], check=True)
    time.sleep(5)

    return Path.cwd().glob("*.fastq")


def handle_available_fastqs(fetch_bundle: FetchBundle) -> None:
    """
    TODO
    """
    output_fasta = f"{fetch_bundle.sra_acc}.collapsed.fa"

    # make sure the file isn't orphaned
    if Path(f"{fetch_bundle.sra_acc}_2.fastq").is_file():
        logger.warning(
            f"""
            An orphaned R2 file was downloaded for {fetch_bundle.sra_acc}, {fetch_bundle.sra_acc}_2.fastq.
            The pipeline's current default behavior is to skip these cases.
            """,  # noqa: E501
        )
        return

    # handle the accession based on whether it contains paired or unpaired reads
    if (
        Path(f"{fetch_bundle.sra_acc}_1.fastq").is_file()
        and Path(f"{fetch_bundle.sra_acc}_2.fastq").is_file()
    ):
        logger.info(f"Processing {fetch_bundle.sra_acc} as paired reads.")

        # merge any overlapping reads, putting unmerged forward and reverse reads into
        # their own FASTQs
        subprocess.run ([
            fetch_bundle.bbmerge_path,
            "qtrim=t",
            f"in1={fetch_bundle.sra_acc}_1.fastq in2={fetch_bundle.sra_acc}_2.fastq",
            f"out={fetch_bundle.sra_acc}.merge.fastq.gz",
            f"outu1={fetch_bundle.sra_acc}.unmerged_r1.fastq.gz",
            f"outu2={fetch_bundle.sra_acc}.unmerged_r2.fastq.gz",
        ], check=True,
        )

        # do away with the original FASTQs
        # Uses path library to handle relative path handling
        Path(f"{fetch_bundle.sra_acc}_1.fastq").unlink(missing_ok=True)
        Path(f"{fetch_bundle.sra_acc}_2.fastq").unlink(missing_ok=True)

        # Combine merged and unmerged reads into a single FASTQ
        subprocess.run(["cat",  # noqa: S603, S607 cat process being called with relative path name
                        f"{fetch_bundle.sra_acc}.merge.fastq.gz",
                        f"{fetch_bundle.sra_acc}.unmerged_r1.fastq.gz",
                        f"{fetch_bundle.sra_acc}.unmerged_r2.fastq.gz",
                        ">",
                        f"{fetch_bundle.sra_acc}.all.fastq.gz",
        ], check=True,
        )

        # clean up clean up everybody everywhere
        Path(f"{fetch_bundle.sra_acc}.merge.fastq.gz").unlink(missing_ok=True)
        Path(f"{fetch_bundle.sra_acc}.unmerged_r1.fastq.gz").unlink(missing_ok=True)
        Path(f"{fetch_bundle.sra_acc}.unmerged_r2.fastq.gz").unlink(missing_ok=True)

        # if a mysterious third FASTQ came along, add that in too
        if Path(f"{fetch_bundle.sra_acc}.fastq").is_file():
            logger.info(
                "Combining all available reads from the three downloaded FASTQs.",
            )
            subprocess.run(["cat",  # noqa: S603, S607
                            f"{fetch_bundle.sra_acc}.fastq",
                            "|",
                            "gzip",
                            "-c",
                            ">>",
                            f"{fetch_bundle.sra_acc}.all.fastq.gz",
                            ], check=True,
                            )
            Path(f"{fetch_bundle.sra_acc}.fastq").unlink(missing_ok=True)

        # Dereplicate the reads
        logger.info("Dereplicating the reads")
        dereplicate(
            f"{fetch_bundle.sra_acc}.all.fastq.gz",
            f"{output_fasta}",
        )
        Path(f"{fetch_bundle.sra_acc}.all.fastq.gz").unlink(missing_ok=True)

    elif (
        Path(f"{fetch_bundle.sra_acc}.fastq").is_file()
        or Path(f"{fetch_bundle.sra_acc}_1.fastq").is_file()
    ):
        logger.info(f"Processing {fetch_bundle.sra_acc} as singleton reads.")
        dereplicate(f"{fetch_bundle.sra_acc}.fastq", f"{output_fasta}")
        Path(f"{fetch_bundle.sra_acc}.fastq").unlink(missing_ok=True)

    elif Path(f"{fetch_bundle.sra_acc}_1.fastq").is_file():
        logger.info(f"Processing {fetch_bundle.sra_acc} as singleton reads.")
        dereplicate(f"{fetch_bundle.sra_acc}_1.fastq", f"{output_fasta}")
        Path(f"{fetch_bundle.sra_acc}.fastq").unlink(missing_ok=True)
    else:
        logger.warning(
            f"Unsupported FASTQ name downloaded for {fetch_bundle.sra_acc}. Skipping",
        )
        return


def align_derep_reads(input_fasta: str, fetch_bundle: FetchBundle) -> None:
    """
    TODO
    """
    # run final mapping to the reference in preparation for haplotype calling
    logger.info(f"Mapping uncompressed file for {fetch_bundle.sra_acc}.")
    subprocess.run(["minimap",  # noqa: S603, S607
                    "-a",
                    f"{fetch_bundle.ref}",
                    f"{input_fasta}",
                    "--sam-hit-only",
                    "--secondary=no",
                    "-o"
                    f"{fetch_bundle.sra_acc}.wg.sam",
                    ], check=True,
                    )

    # Remove the input FASTA
    Path(input_fasta).unlink(missing_ok=True)
    logger.info(f"Finished fetching and deplicating {fetch_bundle.sra_acc}.")
