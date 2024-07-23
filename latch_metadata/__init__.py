
from latch.types.directory import LatchDir
from latch.types.metadata import LatchAuthor, NextflowMetadata, NextflowRuntimeResources

from .parameters import generated_parameters

NextflowMetadata(
    name="Methylseq",
    display_name='nf-core/methylseq',
    author=LatchAuthor(
        name="nf-core",
    ),
    parameters=generated_parameters,
    runtime_resources=NextflowRuntimeResources(
        cpus=4,
        memory=8,
        storage_gib=100,
    ),
    log_dir=LatchDir("latch:///your_log_dir"),
)
