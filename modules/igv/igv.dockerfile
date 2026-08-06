# Auto-generated Dockerfile for hisat2 conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/igv/igv.yaml
# Build: docker build -t hisat2 -f Dockerfile .
# Save:  docker save hisat2 -o /home/luosg/Database/env/hisat2.tar

FROM continuumio/miniconda3:latest

COPY igv.yaml /tmp/igv.yaml

RUN conda env create -f /tmp/igv.yaml && \
    conda clean -afy && \
    rm /tmp/igv.yaml

ENV CONDA_DEFAULT_ENV=hisat2
ENV PATH="/opt/conda/envs/hisat2/bin:$PATH"

CMD ["bash"]
