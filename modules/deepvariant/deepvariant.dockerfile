# Auto-generated Dockerfile for deepvariant conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/deepvariant/deepvariant.yaml
# Build: docker build -t deepvariant -f Dockerfile .
# Save:  docker save deepvariant -o /home/luosg/Database/env/deepvariant.tar

FROM continuumio/miniconda3:latest

COPY deepvariant.yaml /tmp/deepvariant.yaml

RUN conda env create -f /tmp/deepvariant.yaml && \
    conda clean -afy && \
    rm /tmp/deepvariant.yaml

ENV CONDA_DEFAULT_ENV=deepvariant
ENV PATH="/opt/conda/envs/deepvariant/bin:$PATH"

CMD ["bash"]
