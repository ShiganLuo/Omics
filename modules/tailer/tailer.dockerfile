# Auto-generated Dockerfile for tailer conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/tailer/tailer.yaml
# Build: docker build -t tailer -f Dockerfile .
# Save:  docker save tailer -o /home/luosg/Database/env/tailer.tar

FROM continuumio/miniconda3:latest

COPY tailer.yaml /tmp/tailer.yaml

RUN conda env create -f /tmp/tailer.yaml && \
    conda clean -afy && \
    rm /tmp/tailer.yaml

ENV CONDA_DEFAULT_ENV=tailer
ENV PATH="/opt/conda/envs/tailer/bin:$PATH"

CMD ["bash"]
