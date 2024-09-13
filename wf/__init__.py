from typing import Annotated, List, Optional

from flytekit.core.annotation import FlyteAnnotation
from latch.resources.launch_plan import LaunchPlan
from latch.resources.workflow import workflow
from latch.types import metadata
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile

from wf.entrypoint import (
    Aligner,
    Genome,
    Reference_Type,
    Sample,
    initialize,
    nextflow_runtime,
)


@workflow(metadata._nextflow_metadata)
def nf_nf_core_methylseq(
    input: List[Sample],
    run_name: Annotated[
        str,
        FlyteAnnotation(
            {
                "rules": [
                    {
                        "regex": r"^[a-zA-Z0-9_-]+$",
                        "message": "ID name must contain only letters, digits, underscores, and dashes. No spaces are allowed.",
                    }
                ],
            }
        ),
    ],
    genome_source: str,
    genome: Optional[Genome],
    email: Optional[str],
    multiqc_title: Optional[str],
    fasta: Optional[LatchFile],
    fasta_index: Optional[LatchFile],
    bismark_index: Optional[LatchDir],
    bwa_meth_index: Optional[LatchDir],
    # Alignment
    comprehensive: bool,
    non_directional: bool,
    cytosine_report: bool,
    # Special library types
    pbat: bool,
    rrbs: bool,
    slamseq: bool,
    em_seq: bool,
    single_cell: bool,
    accel: bool,
    cegx: bool,
    epignome: bool,
    zymo: bool,
    # Save intermediate files
    save_reference: bool,
    save_align_intermeds: bool,
    unmapped: bool,
    save_trimmed: bool,
    relax_mismatches: bool,
    known_splices: Optional[LatchFile],
    local_alignment: bool,
    minins: Optional[int],
    maxins: Optional[int],
    nomeseq: bool,
    # bwa-meth options
    min_depth: Optional[int],
    ignore_flags: bool,
    methyl_kit: bool,
    # Qualimap Options
    bamqc_regions_file: Optional[LatchFile],
    # Skip pipeline steps
    skip_trimming: bool,
    skip_deduplication: bool,
    skip_multiqc: bool,
    # Additional option
    multiqc_methods_description: Optional[str],
    # Adapter Trimming
    clip_r1: int = 0,
    clip_r2: int = 0,
    three_prime_clip_r1: int = 0,
    three_prime_clip_r2: int = 0,
    nextseq_trim: int = 0,
    # Bismark options
    num_mismatches: float = 0.6,
    meth_cutoff: Optional[int] = None,
    no_overlap: bool = True,
    ignore_r1: Optional[int] = None,
    ignore_r2: int = 2,
    ignore_3prime_r1: Optional[int] = None,
    ignore_3prime_r2: int = 2,
    # Reference Genome
    latch_genome: Reference_Type = Reference_Type.homo_sapiens,
    aligner: Aligner = Aligner.bismark,
    outdir: LatchOutputDir = LatchOutputDir("latch:///Methylseq"),
) -> None:
    """
    nf-core/methylseq is a bioinformatics analysis pipeline used for Methylation (Bisulfite) sequencing data. It pre-processes raw data from FastQ inputs, aligns the reads and performs extensive quality-control on the results.

    <html>
    <p align="center">
    <img src="https://user-images.githubusercontent.com/31255434/182289305-4cc620e3-86ae-480f-9b61-6ca83283caa5.jpg" alt="Latch Verified" width="100">
    </p>

    <p align="center">
    <strong>
    Latch Verified
    </strong>
    </p>

    <p align="center">

    [![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.1343417-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.1343417)

    **nf-core/methylseq** is a bioinformatics analysis pipeline used for Methylation (Bisulfite) sequencing data. It pre-processes raw data from FastQ inputs, aligns the reads and performs extensive quality-control on the results.

    The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker / Singularity containers making installation trivial and results highly reproducible.

    This workflow is hosted on Latch Workflows, using a native Nextflow integration, with a graphical interface for accessible analysis by scientists. There is also an integration with Latch Registry so that batched workflows can be launched from “graphical sample sheets” or tables associating raw sequencing files with metadata.

    ## Pipeline Summary

    The pipeline allows you to choose between running either [Bismark](https://github.com/FelixKrueger/Bismark) or [bwa-meth](https://github.com/brentp/bwa-meth) / [MethylDackel](https://github.com/dpryan79/methyldackel).
    Choose between workflows by using `--aligner bismark` (default, uses bowtie2 for alignment), `--aligner bismark_hisat` or `--aligner bwameth`.

    | Step                                         | Bismark workflow | bwa-meth workflow     |
    | -------------------------------------------- | ---------------- | --------------------- |
    | Generate Reference Genome Index _(optional)_ | Bismark          | bwa-meth              |
    | Merge re-sequenced FastQ files               | cat              | cat                   |
    | Raw data QC                                  | FastQC           | FastQC                |
    | Adapter sequence trimming                    | Trim Galore!     | Trim Galore!          |
    | Align Reads                                  | Bismark          | bwa-meth              |
    | Deduplicate Alignments                       | Bismark          | Picard MarkDuplicates |
    | Extract methylation calls                    | Bismark          | MethylDackel          |
    | Sample report                                | Bismark          | -                     |
    | Summary Report                               | Bismark          | -                     |
    | Alignment QC                                 | Qualimap         | Qualimap              |
    | Sample complexity                            | Preseq           | Preseq                |
    | Project Report                               | MultiQC          | MultiQC               |

    ## Usage

    First, prepare a samplesheet with your input data that looks as follows:

    `samplesheet.csv`:

    sample,fastq_1,fastq_2
    SRR389222_sub1,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub1.fastq.gz
    SRR389222_sub2,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub2.fastq.gz
    SRR389222_sub2,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub3.fastq.gz
    Ecoli_10K_methylated,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/Ecoli_10K_methylated_R1.fastq.gz,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/Ecoli_10K_methylated_R2.fastq.gz


    Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

    ## Pipeline output

    To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/methylseq/results) tab on the nf-core website pipeline page.
    For more details about the output files and reports, please refer to the
    [output documentation](https://nf-co.re/methylseq/output).

    ## Credits

    These scripts were originally written for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/) at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden.

    - Main author:
    - Phil Ewels ([@ewels](https://github.com/ewels/))
    - Maintainers:
    - Felix Krueger ([@FelixKrueger](https://github.com/FelixKrueger))
    - Sateesh Peri ([@Sateesh_Peri](https://github.com/sateeshperi))
    - Edmund Miller ([@EMiller88](https://github.com/emiller88))
    - Contributors:
    - Rickard Hammarén ([@Hammarn](https://github.com/Hammarn/))
    - Alexander Peltzer ([@apeltzer](https://github.com/apeltzer/))
    - Patrick Hüther ([@phue](https://github.com/phue/))

    ## Contributions and Support

    If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

    For further information or help, don't hesitate to get in touch on the [Slack `#methylseq` channel](https://nfcore.slack.com/channels/methylseq) (you can join with [this invite](https://nf-co.re/join/slack)).

    ## Citations

    If you use nf-core/methylseq for your analysis, please cite it using the following doi: [10.5281/zenodo.1343417](https://doi.org/10.5281/zenodo.1343417)

    An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

    You can cite the `nf-core` publication as follows:

    > **The nf-core framework for community-curated bioinformatics pipelines.**
    >
    > Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
    >
    > _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

    """

    pvc_name: str = initialize(run_name=run_name)
    nextflow_runtime(
        pvc_name=pvc_name,
        run_name=run_name,
        input=input,
        outdir=outdir,
        email=email,
        multiqc_title=multiqc_title,
        genome_source=genome_source,
        latch_genome=latch_genome,
        genome=genome,
        fasta=fasta,
        fasta_index=fasta_index,
        bismark_index=bismark_index,
        bwa_meth_index=bwa_meth_index,
        aligner=aligner,
        comprehensive=comprehensive,
        non_directional=non_directional,
        cytosine_report=cytosine_report,
        pbat=pbat,
        rrbs=rrbs,
        slamseq=slamseq,
        em_seq=em_seq,
        single_cell=single_cell,
        accel=accel,
        cegx=cegx,
        epignome=epignome,
        zymo=zymo,
        save_reference=save_reference,
        save_align_intermeds=save_align_intermeds,
        unmapped=unmapped,
        save_trimmed=save_trimmed,
        clip_r1=clip_r1,
        clip_r2=clip_r2,
        three_prime_clip_r1=three_prime_clip_r1,
        three_prime_clip_r2=three_prime_clip_r2,
        nextseq_trim=nextseq_trim,
        relax_mismatches=relax_mismatches,
        num_mismatches=num_mismatches,
        meth_cutoff=meth_cutoff,
        no_overlap=no_overlap,
        ignore_r1=ignore_r1,
        ignore_r2=ignore_r2,
        ignore_3prime_r1=ignore_3prime_r1,
        ignore_3prime_r2=ignore_3prime_r2,
        known_splices=known_splices,
        local_alignment=local_alignment,
        minins=minins,
        maxins=maxins,
        nomeseq=nomeseq,
        min_depth=min_depth,
        ignore_flags=ignore_flags,
        methyl_kit=methyl_kit,
        bamqc_regions_file=bamqc_regions_file,
        skip_trimming=skip_trimming,
        skip_deduplication=skip_deduplication,
        skip_multiqc=skip_multiqc,
        multiqc_methods_description=multiqc_methods_description,
    )


LaunchPlan(
    nf_nf_core_methylseq,
    "Test Data",
    {
        "input": [
            Sample(
                sample="SRR389222_sub1",
                fastq_1=LatchFile(
                    "s3://latch-public/nf-core/methylseq/test_data/SRR389222_sub1.fastq.gz"
                ),
                fastq_2=None,
            ),
            Sample(
                sample="SRR389222_sub2",
                fastq_1=LatchFile(
                    "s3://latch-public/nf-core/methylseq/test_data/SRR389222_sub2.fastq.gz"
                ),
                fastq_2=None,
            ),
            Sample(
                sample="SRR389222_sub3",
                fastq_1=LatchFile(
                    "s3://latch-public/nf-core/methylseq/test_data/SRR389222_sub3.fastq.gz"
                ),
                fastq_2=None,
            ),
            Sample(
                sample="Ecoli_10K_methylated",
                fastq_1=LatchFile(
                    "s3://latch-public/nf-core/methylseq/test_data/Ecoli_10K_methylated_R1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/nf-core/methylseq/test_data/Ecoli_10K_methylated_R2.fastq.gz"
                ),
            ),
        ],
        "run_name": "Test_Run",
        "genome_source": "custom",
        "fasta": LatchFile("s3://latch-public/nf-core/methylseq/test_data/genome.fa"),
        "fasta_index": LatchFile(
            "s3://latch-public/nf-core/methylseq/test_data/genome.fa.fai"
        ),
    },
)
