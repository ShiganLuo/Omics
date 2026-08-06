# Auto-generated Dockerfile for SOAPnuke conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/SOAPnuke/SOAPnuke.yaml
# Build: docker build -t SOAPnuke -f Dockerfile .
# Save:  docker save SOAPnuke -o /home/luosg/Database/env/SOAPnuke.tar

FROM continuumio/miniconda3:latest

COPY SOAPnuke.yaml /tmp/SOAPnuke.yaml

RUN conda env create -f /tmp/SOAPnuke.yaml && \
    conda clean -afy && \
    rm /tmp/SOAPnuke.yaml

ENV CONDA_DEFAULT_ENV=SOAPnuke
ENV PATH="/opt/conda/envs/SOAPnuke/bin:$PATH"

CMD ["bash"]
