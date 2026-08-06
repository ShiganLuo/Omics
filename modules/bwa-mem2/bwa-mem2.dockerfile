# Auto-generated Dockerfile for bwa-mem2 conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/bwa-mem2/bwa-mem2.yaml
# Build: docker build -t bwa-mem2 -f Dockerfile .
# Save:  docker save bwa-mem2 -o /home/luosg/Database/env/bwa-mem2.tar

FROM continuumio/miniconda3:latest

COPY bwa-mem2.yaml /tmp/bwa-mem2.yaml

RUN conda env create -f /tmp/bwa-mem2.yaml && \
    conda clean -afy && \
    rm /tmp/bwa-mem2.yaml

ENV CONDA_DEFAULT_ENV=bwa-mem2
ENV PATH="/opt/conda/envs/bwa-mem2/bin:$PATH"

CMD ["bash"]
