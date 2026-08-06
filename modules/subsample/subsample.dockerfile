# Auto-generated Dockerfile for subsample conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/subsample/subsample.yaml
# Build: docker build -t subsample -f Dockerfile .
# Save:  docker save subsample -o /home/luosg/Database/env/subsample.tar

FROM continuumio/miniconda3:latest

COPY subsample.yaml /tmp/subsample.yaml

RUN conda env create -f /tmp/subsample.yaml && \
    conda clean -afy && \
    rm /tmp/subsample.yaml

ENV CONDA_DEFAULT_ENV=subsample
ENV PATH="/opt/conda/envs/subsample/bin:$PATH"

CMD ["bash"]
