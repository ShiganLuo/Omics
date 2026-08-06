# Auto-generated Dockerfile for arriba conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/arriba/arriba.yaml
# Build: docker build -t arriba -f Dockerfile .
# Save:  docker save arriba -o /home/luosg/Database/env/arriba.tar

FROM continuumio/miniconda3:latest

COPY arriba.yaml /tmp/arriba.yaml

RUN conda env create -f /tmp/arriba.yaml && \
    conda clean -afy && \
    rm /tmp/arriba.yaml

ENV CONDA_DEFAULT_ENV=arriba
ENV PATH="/opt/conda/envs/arriba/bin:$PATH"

CMD ["bash"]
