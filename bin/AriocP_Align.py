import os
import subprocess
import sys
from pathlib import Path

reads_1 = sys.argv[1]
reads_2 = sys.argv[2]
arioc_idx = sys.argv[3]
report_file = "AriocP.report.txt"

bam_output_file = "AriocP.aligned.bam"


def Encode_Reads(reads_1, reads_2):
    data_dir = str(Path().resolve())
    reads_1_name = Path(reads_1).name.replace(".gz", "")
    reads_2_name = Path(reads_2).name.replace(".gz", "")

    reads_1_uncompressed = f"{data_dir}{reads_1_name}"
    reads_2_uncompressed = f"{data_dir}{reads_2_name}"

    if reads_1.endswith(".gz"):
        subprocess.run(f"gunzip -c {reads_1} > {reads_1_uncompressed}", shell=True)
    if reads_1.endswith(".gz"):
        subprocess.run(f"gunzip -c {reads_2} > {reads_2_uncompressed}", shell=True)

    encode_paired_end_reads = f"""
	<?xml version="1.0" encoding="utf-8"?>
	    <AriocE>
	        <dataIn sequenceType="Q">
	            <file subId="1" mate="1">{reads_1_uncompressed}</file>
	            <file subId="1" mate="2">{reads_2_uncompressed}</file>
	         </dataIn>
	         <dataOut>
	             <path>{data_dir}/encoded_reads</path>
	         </dataOut>
	    </AriocE>"""
    f = open("encode_paired_reads.cfg", "w")
    f.write(encode_paired_end_reads)
    f.close()

    subprocess.run(f"./AriocE encode_paired_reads.cfg", shell=True)
    read_filename = Path(reads_1_uncompressed).name
    os.remove(reads_1_uncompressed)
    os.remove(reads_1_uncompressed)

    return f"{data_dir}/encoded_reads", read_filename


def Align_Reads(
    encoded_reads_path, arioc_idx, bam_output_file, report_file, read_filename
):
    temp_sam = "Alignments.sam"
    alignment_cfg = f"""
    <?xml version="1.0" encoding="utf-8"?>
        <AriocP gpuMask="0x0000000F" batchSize="20k">
             <R>{arioc_idx}</R>
             <nongapped seed="ssi84_2_30_CT" maxJ="*" />
             <gapped seed="hsi25_0_30_CT" maxJ="16" Wmxgs="2,6,5,3" Vt="150"/>
             <Q filePath="{encoded_reads_path}">
                 <paired subId="1">
                     <file>{read_filename}.R1</file>
                     <file>{read_filename}.R2</file>
                 </paired>
             </Q>

             <A overwrite="true" cigarFormat="MIDS">
                 <sam report="mu">{temp_sam}</sam>
             </A>
        </AriocP>"""
    f = open("Align_Paired_End_Reads.cfg", "w")
    f.write(alignment_cfg)
    f.close()

    subprocess.run(
        f"/usr/bin/time -v ../bin/AriocP Align_Paired_End_Reads.cfg > {report_file} 2>&1",
        shell=True,
    )
    subprocess.run(
        f"samtools view --threads 20 -bS -F 4 -F 8 -F 256 {temp_sam} > {bam_output_file}",
        shell=True,
    )


encoded_reads, read_filename = Encode_Reads(reads_1, reads_2)
Align_Reads(encoded_reads, arioc_idx, bam_output_file, report_file.read_filename)
