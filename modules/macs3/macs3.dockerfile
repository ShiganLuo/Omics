# Auto-generated Dockerfile for macs3 conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/macs3/macs3.yaml
# Build: docker build -t macs3 -f Dockerfile .
# Save:  docker save macs3 -o /home/luosg/Database/env/macs3.tar

FROM continuumio/miniconda3:latest

COPY macs3.yaml /tmp/macs3.yaml

RUN conda env create -f /tmp/macs3.yaml && \
    conda clean -afy && \
    rm /tmp/macs3.yaml

ENV CONDA_DEFAULT_ENV=macs3
ENV PATH="/opt/conda/envs/macs3/bin:$PATH"

CMD ["bash"]
