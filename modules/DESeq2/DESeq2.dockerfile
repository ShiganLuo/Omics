# Auto-generated Dockerfile for DESeq2 conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/DESeq2/DESeq2.yaml
# Build: docker build -t DESeq2 -f Dockerfile .
# Save:  docker save DESeq2 -o /home/luosg/Database/env/DESeq2.tar

FROM continuumio/miniconda3:latest

COPY DESeq2.yaml /tmp/DESeq2.yaml

RUN conda env create -f /tmp/DESeq2.yaml && \
    conda clean -afy && \
    rm /tmp/DESeq2.yaml

ENV CONDA_DEFAULT_ENV=DESeq2
ENV PATH="/opt/conda/envs/DESeq2/bin:$PATH"

CMD ["bash"]
