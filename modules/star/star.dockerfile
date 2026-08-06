# Auto-generated Dockerfile for star conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/star/star.yaml
# Build: docker build -t star -f Dockerfile .
# Save:  docker save star -o /home/luosg/Database/env/star.tar

FROM continuumio/miniconda3:latest

COPY star.yaml /tmp/star.yaml

RUN conda env create -f /tmp/star.yaml && \
    conda clean -afy && \
    rm /tmp/star.yaml

ENV CONDA_DEFAULT_ENV=star
ENV PATH="/opt/conda/envs/star/bin:$PATH"

CMD ["bash"]
