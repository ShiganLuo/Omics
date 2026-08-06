# Auto-generated Dockerfile for openms-searchengine conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/openms/searchengine/openms.yaml
# Build: docker build -t openms-searchengine -f Dockerfile .
# Save:  docker save openms-searchengine -o /home/luosg/Database/env/openms-searchengine.tar

FROM continuumio/miniconda3:latest

COPY openms.yaml /tmp/openms.yaml

RUN conda env create -f /tmp/openms.yaml && \
    conda clean -afy && \
    rm /tmp/openms.yaml

ENV CONDA_DEFAULT_ENV=openms-searchengine
ENV PATH="/opt/conda/envs/openms-searchengine/bin:$PATH"

CMD ["bash"]
