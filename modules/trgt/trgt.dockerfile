# Auto-generated Dockerfile for trgt conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/trgt/trgt.yaml
# Build: docker build -t trgt -f Dockerfile .
# Save:  docker save trgt -o /home/luosg/Database/env/trgt.tar

FROM continuumio/miniconda3:latest

COPY trgt.yaml /tmp/trgt.yaml

RUN conda env create -f /tmp/trgt.yaml && \
    conda clean -afy && \
    rm /tmp/trgt.yaml

ENV CONDA_DEFAULT_ENV=trgt
ENV PATH="/opt/conda/envs/trgt/bin:$PATH"

CMD ["bash"]
