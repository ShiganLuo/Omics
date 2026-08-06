# Auto-generated Dockerfile for hiphase conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/hiphase/hiphase.yaml
# Build: docker build -t hiphase -f Dockerfile .
# Save:  docker save hiphase -o /home/luosg/Database/env/hiphase.tar

FROM continuumio/miniconda3:latest

COPY hiphase.yaml /tmp/hiphase.yaml

RUN conda env create -f /tmp/hiphase.yaml && \
    conda clean -afy && \
    rm /tmp/hiphase.yaml

ENV CONDA_DEFAULT_ENV=hiphase
ENV PATH="/opt/conda/envs/hiphase/bin:$PATH"

CMD ["bash"]
