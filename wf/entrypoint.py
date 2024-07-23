from dataclasses import dataclass
from enum import Enum
import os
import subprocess
import requests
import shutil
from pathlib import Path
import typing
import typing_extensions

from latch.resources.workflow import workflow
from latch.resources.tasks import nextflow_runtime_task, custom_task
from latch.types.file import LatchFile
from latch.types.directory import LatchDir, LatchOutputDir
from latch.ldata.path import LPath
from latch_cli.nextflow.workflow import get_flag
from latch_cli.nextflow.utils import _get_execution_name
from latch_cli.utils import urljoins
from latch.types import metadata
from flytekit.core.annotation import FlyteAnnotation

from latch_cli.services.register.utils import import_module_by_path

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


@dataclass
class Sample:
    sample: str
    fastq_1: LatchFile
    fastq_2: typing.Optional[LatchFile]




input_construct_samplesheet = metadata._nextflow_metadata.parameters['input'].samplesheet_constructor


@nextflow_runtime_task(cpu=4, memory=8, storage_gib=100)
def nextflow_runtime(pvc_name: str, input: typing.List[Sample], outdir: typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})], email: typing.Optional[str], multiqc_title: typing.Optional[str], save_reference: typing.Optional[bool], save_align_intermeds: typing.Optional[bool], unmapped: typing.Optional[bool], save_trimmed: typing.Optional[bool], genome: typing.Optional[str], fasta: typing.Optional[LatchFile], fasta_index: typing.Optional[LatchFile], bismark_index: typing.Optional[str], bwa_meth_index: typing.Optional[str], comprehensive: typing.Optional[bool], pbat: typing.Optional[bool], rrbs: typing.Optional[bool], slamseq: typing.Optional[bool], em_seq: typing.Optional[bool], single_cell: typing.Optional[bool], accel: typing.Optional[bool], cegx: typing.Optional[bool], epignome: typing.Optional[bool], zymo: typing.Optional[bool], non_directional: typing.Optional[bool], cytosine_report: typing.Optional[bool], relax_mismatches: typing.Optional[bool], meth_cutoff: typing.Optional[int], known_splices: typing.Optional[LatchFile], local_alignment: typing.Optional[bool], minins: typing.Optional[int], maxins: typing.Optional[int], nomeseq: typing.Optional[bool], ignore_flags: typing.Optional[bool], methyl_kit: typing.Optional[bool], bamqc_regions_file: typing.Optional[LatchFile], skip_trimming: typing.Optional[bool], skip_deduplication: typing.Optional[bool], skip_multiqc: typing.Optional[bool], multiqc_methods_description: typing.Optional[str], aligner: str, clip_r1: typing.Optional[int], clip_r2: typing.Optional[int], three_prime_clip_r1: typing.Optional[int], three_prime_clip_r2: typing.Optional[int], nextseq_trim: typing.Optional[int], num_mismatches: typing.Optional[float], no_overlap: typing.Optional[bool], ignore_r1: typing.Optional[int], ignore_r2: typing.Optional[int], ignore_3prime_r1: typing.Optional[int], ignore_3prime_r2: typing.Optional[int], min_depth: typing.Optional[int]) -> None:
    try:
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

        profile_list = ['docker', 'test']
        if False:
            profile_list.extend([p.value for p in execution_profiles])

        if len(profile_list) == 0:
            profile_list.append("standard")

        profiles = ','.join(profile_list)

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
        *get_flag('input', input_samplesheet),
                *get_flag('outdir', outdir),
                *get_flag('email', email),
                *get_flag('multiqc_title', multiqc_title),
                *get_flag('save_reference', save_reference),
                *get_flag('save_align_intermeds', save_align_intermeds),
                *get_flag('unmapped', unmapped),
                *get_flag('save_trimmed', save_trimmed),
                *get_flag('genome', genome),
                *get_flag('fasta', fasta),
                *get_flag('fasta_index', fasta_index),
                *get_flag('bismark_index', bismark_index),
                *get_flag('bwa_meth_index', bwa_meth_index),
                *get_flag('aligner', aligner),
                *get_flag('comprehensive', comprehensive),
                *get_flag('pbat', pbat),
                *get_flag('rrbs', rrbs),
                *get_flag('slamseq', slamseq),
                *get_flag('em_seq', em_seq),
                *get_flag('single_cell', single_cell),
                *get_flag('accel', accel),
                *get_flag('cegx', cegx),
                *get_flag('epignome', epignome),
                *get_flag('zymo', zymo),
                *get_flag('clip_r1', clip_r1),
                *get_flag('clip_r2', clip_r2),
                *get_flag('three_prime_clip_r1', three_prime_clip_r1),
                *get_flag('three_prime_clip_r2', three_prime_clip_r2),
                *get_flag('nextseq_trim', nextseq_trim),
                *get_flag('non_directional', non_directional),
                *get_flag('cytosine_report', cytosine_report),
                *get_flag('relax_mismatches', relax_mismatches),
                *get_flag('num_mismatches', num_mismatches),
                *get_flag('meth_cutoff', meth_cutoff),
                *get_flag('no_overlap', no_overlap),
                *get_flag('ignore_r1', ignore_r1),
                *get_flag('ignore_r2', ignore_r2),
                *get_flag('ignore_3prime_r1', ignore_3prime_r1),
                *get_flag('ignore_3prime_r2', ignore_3prime_r2),
                *get_flag('known_splices', known_splices),
                *get_flag('local_alignment', local_alignment),
                *get_flag('minins', minins),
                *get_flag('maxins', maxins),
                *get_flag('nomeseq', nomeseq),
                *get_flag('min_depth', min_depth),
                *get_flag('ignore_flags', ignore_flags),
                *get_flag('methyl_kit', methyl_kit),
                *get_flag('bamqc_regions_file', bamqc_regions_file),
                *get_flag('skip_trimming', skip_trimming),
                *get_flag('skip_deduplication', skip_deduplication),
                *get_flag('skip_multiqc', skip_multiqc),
                *get_flag('multiqc_methods_description', multiqc_methods_description)
        ]

        print("Launching Nextflow Runtime")
        print(' '.join(cmd))
        print(flush=True)

        env = {
            **os.environ,
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
    finally:
        print()

        nextflow_log = shared_dir / ".nextflow.log"
        if nextflow_log.exists():
            name = _get_execution_name()
            if name is None:
                print("Skipping logs upload, failed to get execution name")
            else:
                remote = LPath(urljoins("latch:///your_log_dir/Methylseq", name, "nextflow.log"))
                print(f"Uploading .nextflow.log to {remote.path}")
                remote.upload_from(nextflow_log)



@workflow(metadata._nextflow_metadata)
def Methylseq(input: typing.List[Sample], outdir: typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})], email: typing.Optional[str], multiqc_title: typing.Optional[str], save_reference: typing.Optional[bool], save_align_intermeds: typing.Optional[bool], unmapped: typing.Optional[bool], save_trimmed: typing.Optional[bool], genome: typing.Optional[str], fasta: typing.Optional[LatchFile], fasta_index: typing.Optional[LatchFile], bismark_index: typing.Optional[str], bwa_meth_index: typing.Optional[str], comprehensive: typing.Optional[bool], pbat: typing.Optional[bool], rrbs: typing.Optional[bool], slamseq: typing.Optional[bool], em_seq: typing.Optional[bool], single_cell: typing.Optional[bool], accel: typing.Optional[bool], cegx: typing.Optional[bool], epignome: typing.Optional[bool], zymo: typing.Optional[bool], non_directional: typing.Optional[bool], cytosine_report: typing.Optional[bool], relax_mismatches: typing.Optional[bool], meth_cutoff: typing.Optional[int], known_splices: typing.Optional[LatchFile], local_alignment: typing.Optional[bool], minins: typing.Optional[int], maxins: typing.Optional[int], nomeseq: typing.Optional[bool], ignore_flags: typing.Optional[bool], methyl_kit: typing.Optional[bool], bamqc_regions_file: typing.Optional[LatchFile], skip_trimming: typing.Optional[bool], skip_deduplication: typing.Optional[bool], skip_multiqc: typing.Optional[bool], multiqc_methods_description: typing.Optional[str], aligner: str = 'bismark', clip_r1: typing.Optional[int] = 0, clip_r2: typing.Optional[int] = 0, three_prime_clip_r1: typing.Optional[int] = 0, three_prime_clip_r2: typing.Optional[int] = 0, nextseq_trim: typing.Optional[int] = 0, num_mismatches: typing.Optional[float] = 0.6, no_overlap: typing.Optional[bool] = True, ignore_r1: typing.Optional[int] = 0, ignore_r2: typing.Optional[int] = 2, ignore_3prime_r1: typing.Optional[int] = 0, ignore_3prime_r2: typing.Optional[int] = 2, min_depth: typing.Optional[int] = 0) -> None:
    """
    nf-core/methylseq

    Sample Description
    """

    pvc_name: str = initialize()
    nextflow_runtime(pvc_name=pvc_name, input=input, outdir=outdir, email=email, multiqc_title=multiqc_title, save_reference=save_reference, save_align_intermeds=save_align_intermeds, unmapped=unmapped, save_trimmed=save_trimmed, genome=genome, fasta=fasta, fasta_index=fasta_index, bismark_index=bismark_index, bwa_meth_index=bwa_meth_index, aligner=aligner, comprehensive=comprehensive, pbat=pbat, rrbs=rrbs, slamseq=slamseq, em_seq=em_seq, single_cell=single_cell, accel=accel, cegx=cegx, epignome=epignome, zymo=zymo, clip_r1=clip_r1, clip_r2=clip_r2, three_prime_clip_r1=three_prime_clip_r1, three_prime_clip_r2=three_prime_clip_r2, nextseq_trim=nextseq_trim, non_directional=non_directional, cytosine_report=cytosine_report, relax_mismatches=relax_mismatches, num_mismatches=num_mismatches, meth_cutoff=meth_cutoff, no_overlap=no_overlap, ignore_r1=ignore_r1, ignore_r2=ignore_r2, ignore_3prime_r1=ignore_3prime_r1, ignore_3prime_r2=ignore_3prime_r2, known_splices=known_splices, local_alignment=local_alignment, minins=minins, maxins=maxins, nomeseq=nomeseq, min_depth=min_depth, ignore_flags=ignore_flags, methyl_kit=methyl_kit, bamqc_regions_file=bamqc_regions_file, skip_trimming=skip_trimming, skip_deduplication=skip_deduplication, skip_multiqc=skip_multiqc, multiqc_methods_description=multiqc_methods_description)

