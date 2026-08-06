# Auto-generated Dockerfile for samtools conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/samtools/samtools.yaml
# Build: docker build -t samtools -f Dockerfile .
# Save:  docker save samtools -o /home/luosg/Database/env/samtools.tar

FROM continuumio/miniconda3:latest

COPY samtools.yaml /tmp/samtools.yaml

RUN conda env create -f /tmp/samtools.yaml && \
    conda clean -afy && \
    rm /tmp/samtools.yaml

ENV CONDA_DEFAULT_ENV=samtools
ENV PATH="/opt/conda/envs/samtools/bin:$PATH"

CMD ["bash"]
