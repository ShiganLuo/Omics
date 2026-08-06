# Auto-generated Dockerfile for homer conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/homer/homer.yaml
# Build: docker build -t homer -f Dockerfile .
# Save:  docker save homer -o /home/luosg/Database/env/homer.tar

FROM continuumio/miniconda3:latest

COPY homer.yaml /tmp/homer.yaml

RUN conda env create -f /tmp/homer.yaml && \
    conda clean -afy && \
    rm /tmp/homer.yaml

ENV CONDA_DEFAULT_ENV=homer
ENV PATH="/opt/conda/envs/homer/bin:$PATH"

CMD ["bash"]
