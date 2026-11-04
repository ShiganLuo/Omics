# Auto-generated Dockerfile for spectrum conda environment
# Source: /data/16T/luosg/workflow/Omics/modules/spectrum/spectrum.yaml
# Build: docker build -t spectrum -f Dockerfile .
# Save:  docker save spectrum -o /data/16T/Database/env/spectrum.tar

FROM continuumio/miniconda3:latest

COPY spectrum.yaml /tmp/spectrum.yaml

RUN conda env create -f /tmp/spectrum.yaml && \
    conda clean -afy && \
    rm /tmp/spectrum.yaml

ENV CONDA_DEFAULT_ENV=spectrum
ENV PATH="/opt/conda/envs/spectrum/bin:$PATH"

CMD ["bash"]
