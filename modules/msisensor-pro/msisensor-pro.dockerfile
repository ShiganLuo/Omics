# Auto-generated Dockerfile for msisensor-pro conda environment
# Source: modules/msisensor-pro/msisensor-pro.yaml
# Build: docker build -t msisensor-pro -f msisensor-pro.dockerfile .
# Save:  docker save msisensor-pro -o /data/16T/Database/env/msisensor-pro/msisensor-pro.tar

FROM continuumio/miniconda3:latest

COPY msisensor-pro.yaml /tmp/msisensor-pro.yaml

RUN conda env create -f /tmp/msisensor-pro.yaml && \
    conda clean -afy && \
    rm /tmp/msisensor-pro.yaml

ENV CONDA_DEFAULT_ENV=msisensor-pro
ENV PATH="/opt/conda/envs/msisensor-pro/bin:$PATH"

CMD ["bash"]
