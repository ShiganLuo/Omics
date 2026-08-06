# Auto-generated Dockerfile for PureCLIP conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/PureCLIP/PureCLIP.yaml
# Build: docker build -t PureCLIP -f Dockerfile .
# Save:  docker save PureCLIP -o /home/luosg/Database/env/PureCLIP.tar

FROM continuumio/miniconda3:latest

COPY PureCLIP.yaml /tmp/PureCLIP.yaml

RUN conda env create -f /tmp/PureCLIP.yaml && \
    conda clean -afy && \
    rm /tmp/PureCLIP.yaml

ENV CONDA_DEFAULT_ENV=PureCLIP
ENV PATH="/opt/conda/envs/PureCLIP/bin:$PATH"

CMD ["bash"]
