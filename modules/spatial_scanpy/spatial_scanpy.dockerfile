# Auto-generated Dockerfile for spatial_scanpy conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/spatial_scanpy/spatial_scanpy.yaml
# Build: docker build -t spatial_scanpy -f Dockerfile .
# Save:  docker save spatial_scanpy -o /home/luosg/Database/env/spatial_scanpy.tar

FROM continuumio/miniconda3:latest

COPY spatial_scanpy.yaml /tmp/spatial_scanpy.yaml

RUN conda env create -f /tmp/spatial_scanpy.yaml && \
    conda clean -afy && \
    rm /tmp/spatial_scanpy.yaml

ENV CONDA_DEFAULT_ENV=spatial_scanpy
ENV PATH="/opt/conda/envs/spatial_scanpy/bin:$PATH"

CMD ["bash"]
