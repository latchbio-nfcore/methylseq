from dataclasses import dataclass
from enum import Enum
from typing import List, Optional

from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch.types.metadata import (
    Fork,
    ForkBranch,
    LatchRule,
    NextflowParameter,
    Params,
    Section,
    Spoiler,
    Text,
)


@dataclass(frozen=True)
class Sample:
    sample: str
    fastq_1: LatchFile
    fastq_2: Optional[LatchFile]


class Genome(Enum):
    hg38 = "hg38"
    hg19 = "hg19"
    mm10 = "mm10"


class Reference_Type(Enum):
    homo_sapiens = "Homo sapiens (RefSeq GRCh38.p14)"
    mus_musculus = "Mus musculus (RefSeq GRCm39)"
    # rattus_norvegicus = "Rattus norvegicus (RefSeq GRCr8)"


class Aligner(Enum):
    bismark = "bismark"
    # bismark_hisat = "bismark_hisat"
    bwameth = "bwameth"


flow = [
    Section(
        "Input",
        Params(
            "input",
        ),
    ),
    Section(
        "Reference Genome",
        Fork(
            "genome_source",
            "",
            latch_genome_source=ForkBranch(
                "Latch Verified Reference Genome",
                Params(
                    "latch_genome",
                ),
            ),
            custom=ForkBranch(
                "Custom Reference",
                Params(
                    "fasta",
                    "fasta_index",
                    "bismark_index",
                    "bwa_meth_index",
                ),
                Spoiler(
                    "Additional options",
                    Text(
                        "Use iGenomes with caution. The transcriptome and GTF files in iGenomes are vastly out of date with respect to current annotations from Ensembl e.g. human iGenomes annotations are from Ensembl release 75, while the current Ensembl release is 108. Please consider downloading and using a more updated version of your reference genome."
                    ),
                    Params(
                        "genome",
                    ),
                ),
            ),
        ),
    ),
    Section(
        "Aligner",
        Params("aligner", "comprehensive"),
    ),
    Section(
        "Output Directory",
        Params("run_name"),
        Text("Parent directory for outputs"),
        Params("outdir"),
    ),
    Spoiler(
        "Advanced options",
        Spoiler(
            "Adapter Trimming",
            Params(
                "clip_r1",
                "clip_r2",
                "three_prime_clip_r1",
                "three_prime_clip_r2",
                "nextseq_trim",
            ),
        ),
        Spoiler(
            "Special Library Types",
            Params(
                "pbat",
                "rrbs",
                "slamseq",
                "em_seq",
                "single_cell",
                "accel",
                "cegx",
                "epignome",
                "zymo",
            ),
        ),
        Spoiler(
            "Bismark Options",
            Params(
                "non_directional",
                "cytosine_report",
                "relax_mismatches",
                "num_mismatches",
                "meth_cutoff",
                "no_overlap",
                "ignore_r1",
                "ignore_r2",
                "ignore_3prime_r1",
                "ignore_3prime_r2",
                "known_splices",
                "local_alignment",
                "minins",
                "maxins",
                "nomeseq",
            ),
        ),
        Spoiler(
            "BWA-meth Options",
            Params(
                "min_depth",
                "ignore_flags",
                "methyl_kit",
            ),
        ),
        Spoiler(
            "Qualimap Options",
            Params(
                "bamqc_regions_file",
            ),
        ),
        Spoiler(
            "Save Intermediate Files",
            Params(
                "save_reference", "save_align_intermeds", "unmapped", "save_trimmed"
            ),
        ),
        Spoiler(
            "Skip Pipeline Steps",
            Params(
                "skip_trimming",
                "skip_deduplication",
                "skip_multiqc",
            ),
        ),
        Spoiler(
            "MultiQC Options",
            Params(
                "email",
                "multiqc_title",
                "multiqc_methods_description",
            ),
        ),
    ),
]


generated_parameters = {
    "input": NextflowParameter(
        type=List[Sample],
        samplesheet=True,
        samplesheet_type="csv",
        display_name="Samplesheet",
    ),
    "run_name": NextflowParameter(
        type=str,
        display_name="Run Name",
        description="Name of run",
        batch_table_column=True,
        rules=[
            LatchRule(
                regex=r"^[a-zA-Z0-9_-]+$",
                message="Run name must contain only letters, digits, underscores, and dashes. No spaces are allowed.",
            )
        ],
    ),
    "outdir": NextflowParameter(
        type=LatchOutputDir,
        section_title=None,
        display_name="Output Directory",
        description="The output directory where the results will be saved.",
    ),
    "email": NextflowParameter(
        type=Optional[str],
        display_name="Email",
        description="Email address for completion summary.",
    ),
    "multiqc_title": NextflowParameter(
        type=Optional[str],
        display_name="MultiQC Report Title",
        description="MultiQC report title. Printed as page header, used for filename if not otherwise specified.",
    ),
    # Reference Genome
    "genome_source": NextflowParameter(),
    "latch_genome": NextflowParameter(
        type=Reference_Type,
        display_name="Latch Verified Reference Genome",
        description="Name of Latch Verified Reference Genome.",
        default=Reference_Type.homo_sapiens,
    ),
    "genome": NextflowParameter(
        type=Genome,
        default=Genome.hg38,
        display_name="Genome",
        description="iGenome genome reference.",
    ),
    "fasta": NextflowParameter(
        type=Optional[LatchFile],
        default=None,
        display_name="FASTA Reference",
        description="FASTA genome file (Only .fa, .fa.gz, .fasta or .fasta.gz accepted)",
    ),
    "fasta_index": NextflowParameter(
        type=Optional[LatchFile],
        display_name="FASTA Index",
        description="Path to Fasta index file.",
    ),
    "bismark_index": NextflowParameter(
        type=Optional[LatchDir],
        default=None,
        display_name="Bismark Index",
        description="Directory containing a Bismark reference index.",
    ),
    "bwa_meth_index": NextflowParameter(
        type=Optional[str],
        display_name="BWA-meth Index",
        description="bwameth index filename base",
    ),
    # Alignment
    "aligner": NextflowParameter(
        type=Aligner,
        default=Aligner.bismark,
        display_name="Aligner",
        description="Alignment tool to use.",
    ),
    "comprehensive": NextflowParameter(
        type=bool,
        default=False,
        display_name="Comprehensive",
        description="Output information for all cytosine contexts.",
    ),
    "non_directional": NextflowParameter(
        type=bool,
        default=False,
        display_name="Non Directional",
        description="Run alignment against all four possible strands.",
    ),
    "cytosine_report": NextflowParameter(
        type=bool,
        default=False,
        display_name="Cytosine Report",
        description="Output stranded cytosine report, following Bismark's bismark_methylation_extractor step.",
    ),
    # Special library types
    "pbat": NextflowParameter(
        type=bool,
        default=False,
        display_name="PBAT",
        description="Preset for working with PBAT libraries.",
    ),
    "rrbs": NextflowParameter(
        type=bool,
        default=False,
        display_name="RRBS",
        description="Turn on if dealing with MspI digested material.",
    ),
    "slamseq": NextflowParameter(
        type=bool,
        default=False,
        display_name="SLAM-seq",
        description="Run bismark in SLAM-seq mode.",
    ),
    "em_seq": NextflowParameter(
        type=bool,
        default=False,
        display_name="EM-seq",
        description="Preset for EM-seq libraries.",
    ),
    "single_cell": NextflowParameter(
        type=bool,
        default=False,
        display_name="Single-cell",
        description="Trimming preset for single-cell bisulfite libraries.",
    ),
    "accel": NextflowParameter(
        type=bool,
        default=False,
        display_name="Accel",
        description="Trimming preset for the Accel kit.",
    ),
    "cegx": NextflowParameter(
        type=bool,
        default=False,
        display_name="CEGX",
        description="Trimming preset for the CEGX bisulfite kit.",
    ),
    "epignome": NextflowParameter(
        type=bool,
        default=False,
        display_name="Epignome",
        description="Trimming preset for the Epignome kit.",
    ),
    "zymo": NextflowParameter(
        type=bool,
        default=False,
        display_name="Zymo",
        description="Trimming preset for the Zymo kit.",
    ),
    # Save intermediate files
    "save_reference": NextflowParameter(
        type=bool,
        default=False,
        display_name="Save Reference",
        description="Save reference(s) to results directory",
    ),
    "save_align_intermeds": NextflowParameter(
        type=bool,
        default=False,
        display_name="Save Align Intermeds",
        description="Save aligned intermediates to results directory",
    ),
    "unmapped": NextflowParameter(
        type=bool,
        default=False,
        display_name="Unmapped",
        description="Save unmapped reads to FastQ files",
    ),
    "save_trimmed": NextflowParameter(
        type=bool,
        default=False,
        display_name="Save Trimmed",
        description="Save trimmed reads to results directory.",
    ),
    # Adapter Trimming
    "clip_r1": NextflowParameter(
        type=int,
        default=10,
        display_name="Clip R1",
        description="Trim bases from the 5' end of read 1 (or single-end reads).",
    ),
    "clip_r2": NextflowParameter(
        type=int,
        default=10,
        display_name="Clip R2",
        description="Trim bases from the 5' end of read 2 (paired-end only).",
    ),
    "three_prime_clip_r1": NextflowParameter(
        type=int,
        default=10,
        display_name="3' Clip R1",
        description="Trim bases from the 3' end of read 1 AFTER adapter/quality trimming.",
    ),
    "three_prime_clip_r2": NextflowParameter(
        type=int,
        default=10,
        display_name="3' Clip R2",
        description="Trim bases from the 3' end of read 2 AFTER adapter/quality trimming",
    ),
    "nextseq_trim": NextflowParameter(
        type=int,
        default=0,
        display_name="Nextseq Trim",
        description="Trim bases below this quality value from the 3' end of the read, ignoring high-quality G bases",
    ),
    # Bismark options
    "relax_mismatches": NextflowParameter(
        type=bool,
        default=False,
        display_name="Relax Mismatches",
        description="Turn on to relax stringency for alignment (set allowed penalty with --num_mismatches).",
    ),
    "num_mismatches": NextflowParameter(
        type=float,
        default=0.6,
        display_name="Number of Mismatches",
        description="0.6 will allow a penalty of bp * -0.6 - for 100bp reads (bismark default is 0.2)",
    ),
    "meth_cutoff": NextflowParameter(
        type=int,
        display_name="Methylation Cutoff",
        description="Specify a minimum read coverage to report a methylation call",
    ),
    "no_overlap": NextflowParameter(
        type=bool,
        default=True,
        display_name="No Overlap",
        description="Ignore read 2 methylation when it overlaps read 1",
    ),
    "ignore_r1": NextflowParameter(
        type=int,
        display_name="Ignore R1",
        description="Ignore methylation in first n bases of 5' end of R1",
    ),
    "ignore_r2": NextflowParameter(
        type=int,
        default=2,
        display_name="Ignore R2",
        description="Ignore methylation in first n bases of 5' end of R2",
    ),
    "ignore_3prime_r1": NextflowParameter(
        type=int,
        display_name="Ignore 3' R1",
        description="Ignore methylation in last n bases of 3' end of R1",
    ),
    "ignore_3prime_r2": NextflowParameter(
        type=int,
        default=2,
        display_name="Ignore 3' R2",
        description="Ignore methylation in last n bases of 3' end of R2",
    ),
    "known_splices": NextflowParameter(
        type=Optional[LatchFile],
        display_name="Known Splices",
        description="Supply a .gtf file containing known splice sites (bismark_hisat only).",
    ),
    "local_alignment": NextflowParameter(
        type=bool,
        default=False,
        display_name="Local Alignment",
        description="Allow soft-clipping of reads (potentially useful for single-cell experiments).",
    ),
    "minins": NextflowParameter(
        type=int,
        display_name="Min Insert Size",
        description="The minimum insert size for valid paired-end alignments.",
    ),
    "maxins": NextflowParameter(
        type=int,
        display_name="Max Insert Size",
        description="The maximum insert size for valid paired-end alignments.",
    ),
    "nomeseq": NextflowParameter(
        type=bool,
        default=False,
        display_name="NOMe-seq",
        description="Sample is NOMe-seq or NMT-seq. Runs coverage2cytosine.",
    ),
    # bwa-meth options
    "min_depth": NextflowParameter(
        type=int,
        display_name="Min Depth",
        description="Specify a minimum read coverage for MethylDackel to report a methylation call.",
    ),
    "ignore_flags": NextflowParameter(
        type=bool,
        default=False,
        display_name="Ignore Flags",
        description="MethylDackel - ignore SAM flags",
    ),
    "methyl_kit": NextflowParameter(
        type=bool,
        default=False,
        display_name="MethylKit",
        description="Save files for use with methylKit",
    ),
    # Qualimap Options
    "bamqc_regions_file": NextflowParameter(
        type=Optional[LatchFile],
        display_name="BamQC Regions File",
        description="A GFF or BED file containing the target regions which will be passed to Qualimap/Bamqc.",
    ),
    # Skip pipeline steps
    "skip_trimming": NextflowParameter(
        type=bool,
        default=False,
        display_name="Skip Trimming",
        description="Skip read trimming.",
    ),
    "skip_deduplication": NextflowParameter(
        type=bool,
        default=False,
        display_name="Skip Deduplication",
        description="Skip deduplication step after alignment.",
    ),
    "skip_multiqc": NextflowParameter(
        type=bool,
        default=False,
        display_name="Skip MultiQC",
        description="Skip MultiQC",
    ),
    # Additional option
    "multiqc_methods_description": NextflowParameter(
        type=Optional[str],
        display_name="MultiQC Methods Description",
        description="Custom MultiQC yaml file containing HTML including a methods description.",
    ),
}
