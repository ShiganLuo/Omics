# Auto-generated Dockerfile for sigtest conda environment
# Source: /data/16T/luosg/workflow/Omics/modules/sigtest/sigtest.yaml
# Build: docker build -t sigtest -f Dockerfile .
# Save:  docker save sigtest -o /data/16T/Database/env/sigtest.tar

FROM continuumio/miniconda3:latest

COPY sigtest.yaml /tmp/sigtest.yaml

RUN conda env create -f /tmp/sigtest.yaml && \
    conda clean -afy && \
    rm /tmp/sigtest.yaml

ENV CONDA_DEFAULT_ENV=sigtest
ENV PATH="/opt/conda/envs/sigtest/bin:$PATH"

CMD ["bash"]
