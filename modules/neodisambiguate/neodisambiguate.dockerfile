# Auto-generated Dockerfile for neodisambiguate conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/neodisambiguate/neodisambiguate.yaml
# Build: docker build -t neodisambiguate -f Dockerfile .
# Save:  docker save neodisambiguate -o /home/luosg/Database/env/neodisambiguate.tar

FROM continuumio/miniconda3:latest

COPY neodisambiguate.yaml /tmp/neodisambiguate.yaml

RUN conda env create -f /tmp/neodisambiguate.yaml && \
    conda clean -afy && \
    rm /tmp/neodisambiguate.yaml

ENV CONDA_DEFAULT_ENV=neodisambiguate
ENV PATH="/opt/conda/envs/neodisambiguate/bin:$PATH"

CMD ["bash"]
