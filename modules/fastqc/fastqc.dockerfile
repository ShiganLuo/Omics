# Auto-generated Dockerfile for fastqc conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/fastqc/fastqc.yaml
# Build: docker build -t fastqc -f Dockerfile .
# Save:  docker save fastqc -o /home/luosg/Database/env/fastqc.tar

FROM continuumio/miniconda3:latest

COPY fastqc.yaml /tmp/fastqc.yaml

RUN conda env create -f /tmp/fastqc.yaml && \
    conda clean -afy && \
    rm /tmp/fastqc.yaml

ENV CONDA_DEFAULT_ENV=fastqc
ENV PATH="/opt/conda/envs/fastqc/bin:$PATH"

CMD ["bash"]
