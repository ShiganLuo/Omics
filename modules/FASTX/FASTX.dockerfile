# Auto-generated Dockerfile for FASTX conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/FASTX/FASTX.yaml
# Build: docker build -t FASTX -f Dockerfile .
# Save:  docker save FASTX -o /home/luosg/Database/env/FASTX.tar

FROM continuumio/miniconda3:latest

COPY FASTX.yaml /tmp/FASTX.yaml

RUN conda env create -f /tmp/FASTX.yaml && \
    conda clean -afy && \
    rm /tmp/FASTX.yaml

ENV CONDA_DEFAULT_ENV=FASTX
ENV PATH="/opt/conda/envs/FASTX/bin:$PATH"

CMD ["bash"]
