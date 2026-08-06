# Auto-generated Dockerfile for bedtools conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/bedtools/bedtools.yaml
# Build: docker build -t bedtools -f Dockerfile .
# Save:  docker save bedtools -o /home/luosg/Database/env/bedtools.tar

FROM continuumio/miniconda3:latest

COPY bedtools.yaml /tmp/bedtools.yaml

RUN conda env create -f /tmp/bedtools.yaml && \
    conda clean -afy && \
    rm /tmp/bedtools.yaml

ENV CONDA_DEFAULT_ENV=bedtools
ENV PATH="/opt/conda/envs/bedtools/bin:$PATH"

CMD ["bash"]
