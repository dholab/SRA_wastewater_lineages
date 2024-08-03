"""
TODO
"""

from __future__ import annotations

import gc
import gzip
from pathlib import Path
from typing import TYPE_CHECKING

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from loguru import logger

if TYPE_CHECKING:
    import io


def handle_opening(input_path: str) -> io.TextIOWrapper | gzip.GzipFile:
    """
    Opens a FASTQ file, supporting both uncompressed and gzip-compressed formats.

    Args:
        input_path (str): Path to the input FASTQ file (.fastq, .fq, or their .gz variants).

    Returns:
        io.TextIOWrapper | gzip.GzipFile: File handle for the opened file.

    Raises:
        AssertionError: If the file doesn't exist or isn't a recognized FASTQ format.

    The function checks if the file exists and has a valid FASTQ extension,
    then opens and returns the appropriate file handle based on whether
    the file is gzip-compressed or not.
    """
    assert Path(
        input_path,
    ).is_file(), f"The provided file, {input_path}, does not exist."
    assert (
        "fastq" in input_path or "fq" in input_path
    ), f"The provided input is not an accepted file type: {input_path}"

    if input_path.endswith(".gz"):
        with gzip.open(input_path, "rb") as gz_handle:
            return gz_handle

    with Path(input_path).open("r") as handle:
        return handle


def dereplicate(fastq_path: str, derep_output: str, min_count: int = 1) -> None:
    """
    Dereplicates sequences from a FASTQ file and writes unique sequences to a FASTA file.

    Args:
        fastq_path (str): Path to the input FASTQ file.
        derep_output (str): Path to the output FASTA file.
        min_count (int, optional): Minimum occurrence count to keep a sequence. Defaults to 1.

    Returns:
        None

    The function reads sequences from the input FASTQ file, counts occurrences,
    filters by minimum count, sorts by frequency, and writes unique sequences
    to the output FASTA file with IDs indicating their order and count.
    """
    # open the file
    file_buffer = handle_opening(fastq_path)

    # parse out the records, collect them into a dictionary where each record sequence
    # is alongside its number of occurrences, and only keep records that occur more than
    # the min count
    record_parser = SeqIO.parse(file_buffer, "fastq")
    total_reads: int = 0
    seq_count_dict: dict[str, int] = {}
    for i, record in enumerate(record_parser):
        if record.seq not in seq_count_dict:
            seq_count_dict[record.seq] = 1
        else:
            seq_count_dict[record.seq] += 1
        total_reads = i
    checked_counts = {
        seq: count for seq, count in seq_count_dict.items() if count >= min_count
    }

    # close the file buffer
    file_buffer.close()

    logger.info(
        f"{total_reads} sequences collected.  Dereplicated to {len(checked_counts)} unique sequences. Writing output file",  # noqa: E501
    )

    # sort the reads in descending order by how frequent they are
    seqs_by_count = sorted(
        checked_counts.items(),
        key=lambda x: x[1],
        reverse=True,
    )

    # instantiate new FASTA records and clear the memory held by the old reads
    final_records = [
        SeqRecord(Seq(seq), id=f"{i}-{count}")
        for i, (seq, count) in enumerate(seqs_by_count)
    ]
    del seqs_by_count
    gc.collect()

    # export the sorted, filtered reads in the simpler FASTA format
    with Path(derep_output).open("w", encoding="utf8") as output_handle:
        SeqIO.write(final_records, output_handle, "fasta")
