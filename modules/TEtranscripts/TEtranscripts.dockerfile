# Auto-generated Dockerfile for TEtranscripts conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/TEtranscripts/TEtranscripts.yaml
# Build: docker build -t TEtranscripts -f Dockerfile .
# Save:  docker save TEtranscripts -o /home/luosg/Database/env/TEtranscripts.tar

FROM continuumio/miniconda3:latest

COPY TEtranscripts.yaml /tmp/TEtranscripts.yaml

RUN conda env create -f /tmp/TEtranscripts.yaml && \
    conda clean -afy && \
    rm /tmp/TEtranscripts.yaml

ENV CONDA_DEFAULT_ENV=TEtranscripts
ENV PATH="/opt/conda/envs/TEtranscripts/bin:$PATH"

CMD ["bash"]
