# Auto-generated Dockerfile for Flexbar conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/Flexbar/Flexbar.yaml
# Build: docker build -t Flexbar -f Dockerfile .
# Save:  docker save Flexbar -o /home/luosg/Database/env/Flexbar.tar

FROM continuumio/miniconda3:latest

COPY Flexbar.yaml /tmp/Flexbar.yaml

RUN conda env create -f /tmp/Flexbar.yaml && \
    conda clean -afy && \
    rm /tmp/Flexbar.yaml

ENV CONDA_DEFAULT_ENV=Flexbar
ENV PATH="/opt/conda/envs/Flexbar/bin:$PATH"

CMD ["bash"]
