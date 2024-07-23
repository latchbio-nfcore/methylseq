
import typing
from dataclasses import dataclass

import typing_extensions
from flytekit.core.annotation import FlyteAnnotation
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch.types.metadata import NextflowParameter

# Import these into your `__init__.py` file:
#
# from .parameters import generated_parameters

@dataclass(frozen=True)
class Sample:
    sample: str
    fastq_1: LatchFile
    fastq_2: typing.Optional[LatchFile]

generated_parameters = {
    'input': NextflowParameter(
        type=typing.List[Sample],
        samplesheet=True,
        samplesheet_type='csv',
        section_title='Input/output options',
        description='Path to comma-separated file containing information about the samples in the experiment.',
    ),
    'outdir': NextflowParameter(
        type=typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})],
        default=None,
        section_title=None,
        description='The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.',
    ),
    'email': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Email address for completion summary.',
    ),
    'multiqc_title': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='MultiQC report title. Printed as page header, used for filename if not otherwise specified.',
    ),
    'save_reference': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title='Save intermediate files',
        description='Save reference(s) to results directory',
    ),
    'save_align_intermeds': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Save aligned intermediates to results directory',
    ),
    'unmapped': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Bismark only - Save unmapped reads to FastQ files',
    ),
    'save_trimmed': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Save trimmed reads to results directory.',
    ),
    'genome': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title='Reference genome options',
        description='Name of iGenomes reference.',
    ),
    'fasta': NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        description='Path to FASTA genome file',
    ),
    'fasta_index': NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        description='Path to Fasta index file.',
    ),
    'bismark_index': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Path to a directory containing a Bismark reference index.',
    ),
    'bwa_meth_index': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='bwameth index filename base',
    ),
    'aligner': NextflowParameter(
        type=str,
        default='bismark',
        section_title='Alignment options',
        description='Alignment tool to use.',
    ),
    'comprehensive': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Output information for all cytosine contexts.',
    ),
    'pbat': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title='Special library types',
        description='Preset for working with PBAT libraries.',
    ),
    'rrbs': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Turn on if dealing with MspI digested material.',
    ),
    'slamseq': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Run bismark in SLAM-seq mode.',
    ),
    'em_seq': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Preset for EM-seq libraries.',
    ),
    'single_cell': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Trimming preset for single-cell bisulfite libraries.',
    ),
    'accel': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Trimming preset for the Accel kit.',
    ),
    'cegx': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Trimming preset for the CEGX bisulfite kit.',
    ),
    'epignome': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Trimming preset for the Epignome kit.',
    ),
    'zymo': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Trimming preset for the Zymo kit.',
    ),
    'clip_r1': NextflowParameter(
        type=typing.Optional[int],
        default=0,
        section_title='Adapter Trimming',
        description="Trim bases from the 5' end of read 1 (or single-end reads).",
    ),
    'clip_r2': NextflowParameter(
        type=typing.Optional[int],
        default=0,
        section_title=None,
        description="Trim bases from the 5' end of read 2 (paired-end only).",
    ),
    'three_prime_clip_r1': NextflowParameter(
        type=typing.Optional[int],
        default=0,
        section_title=None,
        description="Trim bases from the 3' end of read 1 AFTER adapter/quality trimming.",
    ),
    'three_prime_clip_r2': NextflowParameter(
        type=typing.Optional[int],
        default=0,
        section_title=None,
        description="Trim bases from the 3' end of read 2 AFTER adapter/quality trimming",
    ),
    'nextseq_trim': NextflowParameter(
        type=typing.Optional[int],
        default=0,
        section_title=None,
        description="Trim bases below this quality value from the 3' end of the read, ignoring high-quality G bases",
    ),
    'non_directional': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title='Bismark options',
        description='Run alignment against all four possible strands.',
    ),
    'cytosine_report': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description="Output stranded cytosine report, following Bismark's bismark_methylation_extractor step.",
    ),
    'relax_mismatches': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Turn on to relax stringency for alignment (set allowed penalty with --num_mismatches).',
    ),
    'num_mismatches': NextflowParameter(
        type=typing.Optional[float],
        default=0.6,
        section_title=None,
        description='0.6 will allow a penalty of bp * -0.6 - for 100bp reads (bismark default is 0.2)',
    ),
    'meth_cutoff': NextflowParameter(
        type=typing.Optional[int],
        default=None,
        section_title=None,
        description='Specify a minimum read coverage to report a methylation call',
    ),
    'no_overlap': NextflowParameter(
        type=typing.Optional[bool],
        default=True,
        section_title=None,
        description='Ignore read 2 methylation when it overlaps read 1',
    ),
    'ignore_r1': NextflowParameter(
        type=typing.Optional[int],
        default=0,
        section_title=None,
        description="Ignore methylation in first n bases of 5' end of R1",
    ),
    'ignore_r2': NextflowParameter(
        type=typing.Optional[int],
        default=2,
        section_title=None,
        description="Ignore methylation in first n bases of 5' end of R2",
    ),
    'ignore_3prime_r1': NextflowParameter(
        type=typing.Optional[int],
        default=0,
        section_title=None,
        description="Ignore methylation in last n bases of 3' end of R1",
    ),
    'ignore_3prime_r2': NextflowParameter(
        type=typing.Optional[int],
        default=2,
        section_title=None,
        description="Ignore methylation in last n bases of 3' end of R2",
    ),
    'known_splices': NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        description='Supply a .gtf file containing known splice sites (bismark_hisat only).',
    ),
    'local_alignment': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Allow soft-clipping of reads (potentially useful for single-cell experiments).',
    ),
    'minins': NextflowParameter(
        type=typing.Optional[int],
        default=None,
        section_title=None,
        description='The minimum insert size for valid paired-end alignments.',
    ),
    'maxins': NextflowParameter(
        type=typing.Optional[int],
        default=None,
        section_title=None,
        description='The maximum insert size for valid paired-end alignments.',
    ),
    'nomeseq': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Sample is NOMe-seq or NMT-seq. Runs coverage2cytosine.',
    ),
    'min_depth': NextflowParameter(
        type=typing.Optional[int],
        default=0,
        section_title='bwa-meth options',
        description='Specify a minimum read coverage for MethylDackel to report a methylation call.',
    ),
    'ignore_flags': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='MethylDackel - ignore SAM flags',
    ),
    'methyl_kit': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Save files for use with methylKit',
    ),
    'bamqc_regions_file': NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title='Qualimap Options',
        description='A GFF or BED file containing the target regions which will be passed to Qualimap/Bamqc.',
    ),
    'skip_trimming': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title='Skip pipeline steps',
        description='Skip read trimming.',
    ),
    'skip_deduplication': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip deduplication step after alignment.',
    ),
    'skip_multiqc': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip MultiQC',
    ),
    'multiqc_methods_description': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title='Generic options',
        description='Custom MultiQC yaml file containing HTML including a methods description.',
    ),
}

