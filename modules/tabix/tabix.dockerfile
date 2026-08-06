# Auto-generated Dockerfile for tabix conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/tabix/tabix.yaml
# Build: docker build -t tabix -f Dockerfile .
# Save:  docker save tabix -o /home/luosg/Database/env/tabix.tar

FROM continuumio/miniconda3:latest

COPY tabix.yaml /tmp/tabix.yaml

RUN conda env create -f /tmp/tabix.yaml && \
    conda clean -afy && \
    rm /tmp/tabix.yaml

ENV CONDA_DEFAULT_ENV=tabix
ENV PATH="/opt/conda/envs/tabix/bin:$PATH"

CMD ["bash"]
