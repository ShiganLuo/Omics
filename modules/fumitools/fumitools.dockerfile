# Auto-generated Dockerfile for fumitools conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/fumitools/fumitools.yaml
# Build: docker build -t fumitools -f Dockerfile .
# Save:  docker save fumitools -o /home/luosg/Database/env/fumitools.tar

FROM continuumio/miniconda3:latest

COPY fumitools.yaml /tmp/fumitools.yaml

RUN conda env create -f /tmp/fumitools.yaml && \
    conda clean -afy && \
    rm /tmp/fumitools.yaml

ENV CONDA_DEFAULT_ENV=fumitools
ENV PATH="/opt/conda/envs/fumitools/bin:$PATH"

CMD ["bash"]
