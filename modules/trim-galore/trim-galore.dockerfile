# Auto-generated Dockerfile for trim-galore conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/trim-galore/trim-galore.yaml
# Build: docker build -t trim-galore -f Dockerfile .
# Save:  docker save trim-galore -o /home/luosg/Database/env/trim-galore.tar

FROM continuumio/miniconda3:latest

COPY trim-galore.yaml /tmp/trim-galore.yaml

RUN conda env create -f /tmp/trim-galore.yaml && \
    conda clean -afy && \
    rm /tmp/trim-galore.yaml

ENV CONDA_DEFAULT_ENV=trim-galore
ENV PATH="/opt/conda/envs/trim-galore/bin:$PATH"

CMD ["bash"]
