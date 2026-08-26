# Auto-generated Dockerfile for scanpy conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/scanpy/scanpy.yaml
# Build: docker build -t scanpy -f Dockerfile .
# Save:  docker save scanpy -o /home/luosg/Database/env/scanpy.tar

FROM continuumio/miniconda3:latest

COPY scanpy.yaml /tmp/scanpy.yaml

RUN conda env create -f /tmp/scanpy.yaml && \
    conda clean -afy && \
    rm /tmp/scanpy.yaml

ENV CONDA_DEFAULT_ENV=scanpy
ENV PATH="/opt/conda/envs/scanpy/bin:$PATH"

CMD ["bash"]
