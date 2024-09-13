import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Annotated, List, Optional

import requests
from flytekit.core.annotation import FlyteAnnotation
from latch.executions import rename_current_execution, report_nextflow_used_storage
from latch.ldata.path import LPath
from latch.resources.tasks import custom_task, nextflow_runtime_task
from latch.types import metadata
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch_cli.nextflow.utils import _get_execution_name
from latch_cli.nextflow.workflow import get_flag
from latch_cli.services.register.utils import import_module_by_path
from latch_cli.utils import urljoins

meta = Path("latch_metadata") / "__init__.py"
import_module_by_path(meta)

sys.stdout.reconfigure(line_buffering=True)

input_construct_samplesheet = metadata._nextflow_metadata.parameters[
    "input"
].samplesheet_constructor


@dataclass(frozen=True)
class Sample:
    sample: Annotated[
        str,
        FlyteAnnotation(
            {
                "rules": [
                    {
                        "regex": r"^[^\s]+$",
                        "message": "Sample name cannot contain spaces.",
                    }
                ]
            }
        ),
    ]
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


@custom_task(cpu=0.25, memory=0.5, storage_gib=1)
def initialize(run_name: str) -> str:
    rename_current_execution(str(run_name))

    token = os.environ.get("FLYTE_INTERNAL_EXECUTION_ID")
    if token is None:
        raise RuntimeError("failed to get execution token")

    headers = {"Authorization": f"Latch-Execution-Token {token}"}

    print("Provisioning shared storage volume... ", end="")
    resp = requests.post(
        "http://nf-dispatcher-service.flyte.svc.cluster.local/provision-storage",
        headers=headers,
        json={
            "storage_expiration_hours": 0,
            "version": 2,
        },
    )
    resp.raise_for_status()
    print("Done.")

    return resp.json()["name"]


@nextflow_runtime_task(cpu=4, memory=8, storage_gib=100)
def nextflow_runtime(
    pvc_name: str,
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
    input: List[Sample],
    outdir: LatchOutputDir,
    email: Optional[str],
    multiqc_title: Optional[str],
    # Reference Genome
    genome_source: str,
    latch_genome: Reference_Type,
    genome: Optional[Genome],
    fasta: Optional[LatchFile],
    fasta_index: Optional[LatchFile],
    bismark_index: Optional[LatchDir],
    bwa_meth_index: Optional[LatchDir],
    # Alignment
    aligner: Aligner,
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
    # Adapter Trimming
    clip_r1: int,
    clip_r2: int,
    three_prime_clip_r1: int,
    three_prime_clip_r2: int,
    nextseq_trim: int,
    # Bismark options
    relax_mismatches: bool,
    num_mismatches: float,
    meth_cutoff: Optional[int],
    no_overlap: bool,
    ignore_r1: Optional[int],
    ignore_r2: Optional[int],
    ignore_3prime_r1: Optional[int],
    ignore_3prime_r2: Optional[int],
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
) -> None:
    shared_dir = Path("/nf-workdir")

    input_samplesheet = input_construct_samplesheet(input)

    ignore_list = [
        "latch",
        ".latch",
        ".git",
        "nextflow",
        ".nextflow",
        "work",
        "results",
        "miniconda",
        "anaconda3",
        "mambaforge",
    ]

    shutil.copytree(
        Path("/root"),
        shared_dir,
        ignore=lambda src, names: ignore_list,
        ignore_dangling_symlinks=True,
        dirs_exist_ok=True,
    )

    cmd = [
        "/root/nextflow",
        "run",
        str(shared_dir / "main.nf"),
        "-work-dir",
        str(shared_dir),
        "-profile",
        "docker",
        "-c",
        "latch.config",
        "-resume",
        *get_flag("input", input_samplesheet),
        *get_flag("outdir", LatchOutputDir(f"{outdir.remote_path}/{run_name}")),
        *get_flag("email", email),
        *get_flag("multiqc_title", multiqc_title),
        # Reference Genome
        *get_flag("genome", genome),
        *get_flag("fasta", fasta),
        *get_flag("fasta_index", fasta_index),
        *get_flag("bismark_index", bismark_index),
        *get_flag("bwa_meth_index", bwa_meth_index),
        # Alignment
        *get_flag("aligner", aligner),
        *get_flag("comprehensive", comprehensive),
        *get_flag("non_directional", non_directional),
        *get_flag("cytosine_report", cytosine_report),
        # Special library types
        *get_flag("pbat", pbat),
        *get_flag("rrbs", rrbs),
        *get_flag("slamseq", slamseq),
        *get_flag("em_seq", em_seq),
        *get_flag("single_cell", single_cell),
        *get_flag("accel", accel),
        *get_flag("cegx", cegx),
        *get_flag("epignome", epignome),
        *get_flag("zymo", zymo),
        # Save intermediate files
        *get_flag("save_reference", save_reference),
        *get_flag("save_align_intermeds", save_align_intermeds),
        *get_flag("unmapped", unmapped),
        *get_flag("save_trimmed", save_trimmed),
        # Adapter Trimming
        *get_flag("clip_r1", clip_r1),
        *get_flag("clip_r2", clip_r2),
        *get_flag("three_prime_clip_r1", three_prime_clip_r1),
        *get_flag("three_prime_clip_r2", three_prime_clip_r2),
        *get_flag("nextseq_trim", nextseq_trim),
        # Bismark options
        *get_flag("relax_mismatches", relax_mismatches),
        *get_flag("num_mismatches", num_mismatches),
        *get_flag("meth_cutoff", meth_cutoff),
        *get_flag("no_overlap", no_overlap),
        *get_flag("ignore_r1", ignore_r1),
        *get_flag("ignore_r2", ignore_r2),
        *get_flag("ignore_3prime_r1", ignore_3prime_r1),
        *get_flag("ignore_3prime_r2", ignore_3prime_r2),
        *get_flag("known_splices", known_splices),
        *get_flag("local_alignment", local_alignment),
        *get_flag("minins", minins),
        *get_flag("maxins", maxins),
        *get_flag("nomeseq", nomeseq),
        # bwa-meth options
        *get_flag("min_depth", min_depth),
        *get_flag("ignore_flags", ignore_flags),
        *get_flag("methyl_kit", methyl_kit),
        # Qualimap Options
        *get_flag("bamqc_regions_file", bamqc_regions_file),
        # Skip pipeline steps
        *get_flag("skip_trimming", skip_trimming),
        *get_flag("skip_deduplication", skip_deduplication),
        *get_flag("skip_multiqc", skip_multiqc),
        # Additional option
        *get_flag("multiqc_methods_description", multiqc_methods_description),
    ]

    if genome_source == "latch_genome_source":
        cmd += [
            "--fasta",
            f"s3://latch-public/nf-core/methylseq/{latch_genome.name}/{latch_genome.name}.genomic.fa",
            "--fasta_index",
            f"s3://latch-public/nf-core/methylseq/{latch_genome.name}/{latch_genome.name}.genomic.fa.fai",
            "--bismark_index",
            f"s3://latch-public/nf-core/methylseq/{latch_genome.name}/BismarkIndex",
            "--bwa_meth_index",
            f"s3://latch-public/nf-core/methylseq/{latch_genome.name}/bwameth",
        ]

    print("Launching Nextflow Runtime")
    print(" ".join(cmd))
    print(flush=True)

    failed = False
    try:
        env = {
            **os.environ,
            "NXF_ANSI_LOG": "false",
            "NXF_HOME": "/root/.nextflow",
            "NXF_OPTS": "-Xms1536M -Xmx6144M -XX:ActiveProcessorCount=4",
            "NXF_DISABLE_CHECK_LATEST": "true",
            "NXF_ENABLE_VIRTUAL_THREADS": "false",
        }
        subprocess.run(
            cmd,
            env=env,
            check=True,
            cwd=str(shared_dir),
        )
    except subprocess.CalledProcessError:
        failed = True
    finally:
        print()

        nextflow_log = shared_dir / ".nextflow.log"
        if nextflow_log.exists():
            name = _get_execution_name()
            if name is None:
                print("Skipping logs upload, failed to get execution name")
            else:
                remote = LPath(
                    urljoins(
                        "latch:///your_log_dir/nf_nf_core_methylseq",
                        name,
                        "nextflow.log",
                    )
                )
                print(f"Uploading .nextflow.log to {remote.path}")
                remote.upload_from(nextflow_log)

        print("Computing size of workdir... ", end="")
        try:
            result = subprocess.run(
                ["du", "-sb", str(shared_dir)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=5 * 60,
            )

            size = int(result.stdout.split()[0])
            report_nextflow_used_storage(size)
            print(f"Done. Workdir size: {size / 1024 / 1024 / 1024: .2f} GiB")
        except subprocess.TimeoutExpired:
            print(
                "Failed to compute storage size: Operation timed out after 5 minutes."
            )
        except subprocess.CalledProcessError as e:
            print(f"Failed to compute storage size: {e.stderr}")
        except Exception as e:
            print(f"Failed to compute storage size: {e}")

    if failed:
        sys.exit(1)
