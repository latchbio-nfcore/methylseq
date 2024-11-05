import os
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO

input_fasta = sys.argv[1]
idx_dir = sys.argv[2]


def split_fasta_by_nucleotides(input_fasta, k, output_prefix="split"):
    nucleotide_count = 0
    file_count = 1
    sequences = []

    output_dir = f"{output_prefix}_files/"
    os.makedirs(output_dir, exist_ok=True)

    for record in SeqIO.parse(input_fasta, "fasta"):
        seq_len = len(record.seq)

        if nucleotide_count + seq_len > k:
            output_file = os.path.join(
                output_dir, f"{output_prefix}_{file_count}.fasta"
            )
            SeqIO.write(sequences, output_file, "fasta")
            print(
                f"Written {len(sequences)} sequences with {nucleotide_count} nucleotides to {output_file}"
            )

            # Start a new batch
            file_count += 1
            nucleotide_count = 0
            sequences = []

        # Add the current sequence to the batch and update the count
        sequences.append(record)
        nucleotide_count += seq_len

    # Write any remaining sequences to a new file
    if sequences:
        output_file = os.path.join(output_dir, f"{output_prefix}_{file_count}.fasta")
        SeqIO.write(sequences, output_file, "fasta")
        print(
            f"Written {len(sequences)} sequences with {nucleotide_count} nucleotides to {output_file}"
        )
    return str(Path(output_dir).resolve())


def Encode_References(ref_dir, prefix, idx_dir):
    files = os.listdir(ref_dir)
    s = ""
    ctr = 1
    for f in files:
        if f.startswith(prefix):
            s += f'\n<file SN="(*)" subId="{ctr}">{ref_dir}f</file>'
    build_cfg = f"""
    <?xml version="1.0" encoding="utf-8"?>
    <AriocE seed="ssi84_2_30_CT" maxJ="200">
        <dataIn sequenceType="R">{s}
        </dataIn>
        <dataOut>
            <path>{idx_dir}</path>
        </dataOut>
    </AriocE>
    """
    f = open("Build_Index.cfg", "w")
    f.write(build_cfg)
    f.close()
    subprocess.run(f"/usr/bin/time -v AriocE Build_Index.cfg", shell=True)


# Example usage:
ref_dir = split_fasta_by_nucleotides(input_fasta, k=100000000)
Encode_References(ref_dir, idx_dir)
