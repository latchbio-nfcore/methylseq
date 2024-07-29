import typing
from dataclasses import dataclass
from enum import Enum

from latch.types.directory import LatchOutputDir
from latch.types.file import LatchFile
from latch.types.metadata import FlowBase, NextflowParameter, Params, Section

# Import these into your `__init__.py` file:
#
# from .parameters import generated_parameters


@dataclass(frozen=True)
class Sample:
    sample: str
    fastq_1: LatchFile
    fastq_2: typing.Optional[LatchFile]


class Genome(Enum):
    hg38 = "hg38"
    hg19 = "hg19"
    mm10 = "mm10"


flow: typing.List[FlowBase] = [
    Section("Input / Output", Params("input", "outdir")),
    Params("genome"),
    Section("Bismark Alignment Options", Params("comprehensive", "non_directional", "cytosine_report")),
    Section(
        "Adapter Trimming", Params("clip_r1", "clip_r2", "three_prime_clip_r1", "three_prime_clip_r2", "nextseq_trim")
    ),
    Section("Save Intermediate Files", Params("save_reference", "save_align_intermeds", "unmapped", "save_trimmed")),
]


generated_parameters = {
    "input": NextflowParameter(
        type=typing.List[Sample], samplesheet=True, samplesheet_type="csv", display_name="Samplesheet"
    ),
    "outdir": NextflowParameter(
        type=LatchOutputDir,
        section_title=None,
        display_name="Output Directory",
        description="The output directory where the results will be saved",
    ),
    # Reference Genome
    # todo(ayush): fork here to allow using custom genomes
    "genome": NextflowParameter(
        type=Genome,
        default=Genome.hg38,
        display_name="Genome",
    ),
    # Alignment
    "comprehensive": NextflowParameter(
        type=bool,
        default=False,
        display_name="Output information for all cytosine contexts.",
    ),
    "non_directional": NextflowParameter(
        type=bool,
        default=False,
        display_name="Run alignment against all four possible strands.",
    ),
    "cytosine_report": NextflowParameter(
        type=bool,
        default=False,
        display_name="Output stranded cytosine report, following Bismark's bismark_methylation_extractor step.",
    ),
    # 'Save intermediate files'
    "save_reference": NextflowParameter(
        type=bool,
        default=False,
        display_name="Save reference(s) to results directory",
    ),
    "save_align_intermeds": NextflowParameter(
        type=bool,
        default=False,
        display_name="Save aligned intermediates to results directory",
    ),
    "unmapped": NextflowParameter(
        type=bool,
        default=False,
        display_name="Save unmapped reads to FastQ files",
    ),
    "save_trimmed": NextflowParameter(
        type=bool,
        default=False,
        display_name="Save trimmed reads to results directory.",
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
}
