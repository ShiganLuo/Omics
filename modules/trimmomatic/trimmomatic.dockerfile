# Auto-generated Dockerfile for trimmomatic conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/trimmomatic/trimmomatic.yaml
# Build: docker build -t trimmomatic -f Dockerfile .
# Save:  docker save trimmomatic -o /home/luosg/Database/env/trimmomatic.tar

FROM continuumio/miniconda3:latest

COPY trimmomatic.yaml /tmp/trimmomatic.yaml

RUN conda env create -f /tmp/trimmomatic.yaml && \
    conda clean -afy && \
    rm /tmp/trimmomatic.yaml

ENV CONDA_DEFAULT_ENV=trimmomatic
ENV PATH="/opt/conda/envs/trimmomatic/bin:$PATH"

CMD ["bash"]
