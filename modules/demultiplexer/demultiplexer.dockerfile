# Auto-generated Dockerfile for demultiplexer conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/demultiplexer/demultiplexer.yaml
# Build: docker build -t demultiplexer -f Dockerfile .
# Save:  docker save demultiplexer -o /home/luosg/Database/env/demultiplexer.tar

FROM continuumio/miniconda3:latest

COPY demultiplexer.yaml /tmp/demultiplexer.yaml

RUN conda env create -f /tmp/demultiplexer.yaml && \
    conda clean -afy && \
    rm /tmp/demultiplexer.yaml

ENV CONDA_DEFAULT_ENV=demultiplexer
ENV PATH="/opt/conda/envs/demultiplexer/bin:$PATH"

CMD ["bash"]
