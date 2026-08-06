# Auto-generated Dockerfile for function conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/function/function.yaml
# Build: docker build -t function -f Dockerfile .
# Save:  docker save function -o /home/luosg/Database/env/function.tar

FROM continuumio/miniconda3:latest

COPY function.yaml /tmp/function.yaml

RUN conda env create -f /tmp/function.yaml && \
    conda clean -afy && \
    rm /tmp/function.yaml

ENV CONDA_DEFAULT_ENV=function
ENV PATH="/opt/conda/envs/function/bin:$PATH"

CMD ["bash"]
