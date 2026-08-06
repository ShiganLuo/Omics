# Auto-generated Dockerfile for RmrRNA conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/RmrRNA/RmrRNA.yaml
# Build: docker build -t RmrRNA -f Dockerfile .
# Save:  docker save RmrRNA -o /home/luosg/Database/env/RmrRNA.tar

FROM continuumio/miniconda3:latest

COPY RmrRNA.yaml /tmp/RmrRNA.yaml

RUN conda env create -f /tmp/RmrRNA.yaml && \
    conda clean -afy && \
    rm /tmp/RmrRNA.yaml

ENV CONDA_DEFAULT_ENV=RmrRNA
ENV PATH="/opt/conda/envs/RmrRNA/bin:$PATH"

CMD ["bash"]
