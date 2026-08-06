# Auto-generated Dockerfile for openms-msstats conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/openms/msstats/openms.yaml
# Build: docker build -t openms-msstats -f Dockerfile .
# Save:  docker save openms-msstats -o /home/luosg/Database/env/openms-msstats.tar

FROM continuumio/miniconda3:latest

COPY openms.yaml /tmp/openms.yaml

RUN conda env create -f /tmp/openms.yaml && \
    conda clean -afy && \
    rm /tmp/openms.yaml

ENV CONDA_DEFAULT_ENV=openms-msstats
ENV PATH="/opt/conda/envs/openms-msstats/bin:$PATH"

CMD ["bash"]
