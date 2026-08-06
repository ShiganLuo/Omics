# Auto-generated Dockerfile for manta conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/manta/manta.yaml
# Build: docker build -t manta -f Dockerfile .
# Save:  docker save manta -o /home/luosg/Database/env/manta.tar

FROM continuumio/miniconda3:latest

COPY manta.yaml /tmp/manta.yaml

RUN conda env create -f /tmp/manta.yaml && \
    conda clean -afy && \
    rm /tmp/manta.yaml

ENV CONDA_DEFAULT_ENV=manta
ENV PATH="/opt/conda/envs/manta/bin:$PATH"

CMD ["bash"]
