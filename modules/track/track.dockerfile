# Auto-generated Dockerfile for track conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/track/track.yaml
# Build: docker build -t track -f Dockerfile .
# Save:  docker save track -o /home/luosg/Database/env/track.tar

FROM continuumio/miniconda3:latest

COPY track.yaml /tmp/track.yaml

RUN conda env create -f /tmp/track.yaml && \
    conda clean -afy && \
    rm /tmp/track.yaml

ENV CONDA_DEFAULT_ENV=track
ENV PATH="/opt/conda/envs/track/bin:$PATH"

CMD ["bash"]
