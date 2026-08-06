# Auto-generated Dockerfile for cutadapt conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/cutadapt/cutadapt.yaml
# Build: docker build -t cutadapt -f Dockerfile .
# Save:  docker save cutadapt -o /home/luosg/Database/env/cutadapt.tar

FROM continuumio/miniconda3:latest

COPY cutadapt.yaml /tmp/cutadapt.yaml

RUN conda env create -f /tmp/cutadapt.yaml && \
    conda clean -afy && \
    rm /tmp/cutadapt.yaml

ENV CONDA_DEFAULT_ENV=cutadapt
ENV PATH="/opt/conda/envs/cutadapt/bin:$PATH"

CMD ["bash"]
