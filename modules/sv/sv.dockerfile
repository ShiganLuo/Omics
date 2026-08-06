# Auto-generated Dockerfile for sv conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/sv/sv.yaml
# Build: docker build -t sv -f Dockerfile .
# Save:  docker save sv -o /home/luosg/Database/env/sv.tar

FROM continuumio/miniconda3:latest

COPY sv.yaml /tmp/sv.yaml

RUN conda env create -f /tmp/sv.yaml && \
    conda clean -afy && \
    rm /tmp/sv.yaml

ENV CONDA_DEFAULT_ENV=sv
ENV PATH="/opt/conda/envs/sv/bin:$PATH"

CMD ["bash"]
