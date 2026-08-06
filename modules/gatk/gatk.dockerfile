# Auto-generated Dockerfile for gatk4 conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/gatk/gatk.yaml
# Build: docker build -t gatk4 -f Dockerfile .
# Save:  docker save gatk4 -o /home/luosg/Database/env/gatk4.tar

FROM continuumio/miniconda3:latest

COPY gatk.yaml /tmp/gatk.yaml

RUN conda env create -f /tmp/gatk.yaml && \
    conda clean -afy && \
    rm /tmp/gatk.yaml

ENV CONDA_DEFAULT_ENV=gatk4
ENV PATH="/opt/conda/envs/gatk4/bin:$PATH"

CMD ["bash"]
