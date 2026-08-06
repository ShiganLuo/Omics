# Auto-generated Dockerfile for openms conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/openms/decoydatabase/openms.yaml
# Build: docker build -t openms -f Dockerfile .
# Save:  docker save openms -o /home/luosg/Database/env/openms.tar

FROM continuumio/miniconda3:latest

COPY openms.yaml /tmp/openms.yaml

RUN conda env create -f /tmp/openms.yaml && \
    conda clean -afy && \
    rm /tmp/openms.yaml

ENV CONDA_DEFAULT_ENV=openms
ENV PATH="/opt/conda/envs/openms/bin:$PATH"

CMD ["bash"]
