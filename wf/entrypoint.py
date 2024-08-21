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
from latch.resources.workflow import workflow
from latch.types import metadata
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch_cli.nextflow.utils import _get_execution_name
from latch_cli.nextflow.workflow import get_flag
from latch_cli.services.register.utils import import_module_by_path
from latch_cli.utils import urljoins

meta = Path("latch_metadata") / "__init__.py"
import_module_by_path(meta)
import latch_metadata


@custom_task(cpu=0.25, memory=0.5, storage_gib=1)
def initialize() -> str:
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
        },
    )
    resp.raise_for_status()
    print("Done.")

    return resp.json()["name"]


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


input_construct_samplesheet = metadata._nextflow_metadata.parameters[
    "input"
].samplesheet_constructor


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
    genome_source: str,
    fasta: Optional[LatchFile],
    bismark_index: Optional[LatchDir],
    genome: Genome,
    comprehensive: bool,
    non_directional: bool,
    cytosine_report: bool,
    save_reference: bool,
    save_align_intermeds: bool,
    unmapped: bool,
    save_trimmed: bool,
    clip_r1: int,
    clip_r2: int,
    three_prime_clip_r1: int,
    three_prime_clip_r2: int,
    nextseq_trim: int,
) -> None:
    shared_dir = Path("/nf-workdir")

    rename_current_execution(str(run_name))

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

    profile_list = []
    # if False:
    #     profile_list.extend([p.value for p in execution_profiles])

    if len(profile_list) == 0:
        profile_list.append("standard")

    profiles = ",".join(profile_list)

    cmd = [
        "/root/nextflow",
        "run",
        str(shared_dir / "main.nf"),
        "-work-dir",
        str(shared_dir),
        "-profile",
        profiles,
        "-c",
        "latch.config",
        "-resume",
        *get_flag("input", input_samplesheet),
        *get_flag("outdir", LatchOutputDir(f"{outdir.remote_path}/{run_name}")),
        *get_flag("comprehensive", comprehensive),
        *get_flag("non_directional", non_directional),
        *get_flag("cytosine_report", cytosine_report),
        *get_flag("save_reference", save_reference),
        *get_flag("save_align_intermeds", save_align_intermeds),
        *get_flag("unmapped", unmapped),
        *get_flag("save_trimmed", save_trimmed),
        *get_flag("clip_r1", clip_r1),
        *get_flag("clip_r2", clip_r2),
        *get_flag("three_prime_clip_r1", three_prime_clip_r1),
        *get_flag("three_prime_clip_r2", three_prime_clip_r2),
        *get_flag("nextseq_trim", nextseq_trim),
    ]

    if genome_source == "custom":
        cmd += [
            *get_flag("fasta", fasta),
            # *get_flag("fasta_index", fasta_index),
            *get_flag("bismark_index", bismark_index),
        ]
    elif genome_source == "igenome":
        cmd += [
            *get_flag("genome", genome),
        ]
    # todo() - compile Latch specific genomes
    # elif genome_source == "latch_genome":

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
                    urljoins("latch:///methylseq-logs/Methylseq", name, "nextflow.log")
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


@workflow(metadata._nextflow_metadata)
def Methylseq(
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
    outdir: LatchOutputDir,
    genome_source: str,
    fasta: Optional[LatchFile] = None,
    bismark_index: Optional[LatchDir] = None,
    genome: Genome = Genome.hg38,
    comprehensive: bool = False,
    non_directional: bool = False,
    cytosine_report: bool = False,
    save_reference: bool = False,
    save_align_intermeds: bool = False,
    unmapped: bool = False,
    save_trimmed: bool = False,
    clip_r1: int = 10,
    clip_r2: int = 10,
    three_prime_clip_r1: int = 10,
    three_prime_clip_r2: int = 10,
    nextseq_trim: int = 0,
) -> None:
    """
    nf-core/methylseq

    Sample Description
    """

    pvc_name: str = initialize()
    nextflow_runtime(
        pvc_name=pvc_name,
        run_name=run_name,
        input=input,
        outdir=outdir,
        genome_source=genome_source,
        fasta=fasta,
        bismark_index=bismark_index,
        genome=genome,
        comprehensive=comprehensive,
        non_directional=non_directional,
        cytosine_report=cytosine_report,
        save_reference=save_reference,
        save_align_intermeds=save_align_intermeds,
        unmapped=unmapped,
        save_trimmed=save_trimmed,
        clip_r1=clip_r1,
        clip_r2=clip_r2,
        three_prime_clip_r1=three_prime_clip_r1,
        three_prime_clip_r2=three_prime_clip_r2,
        nextseq_trim=nextseq_trim,
    )
