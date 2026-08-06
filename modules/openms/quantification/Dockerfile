# Auto-generated Dockerfile for openms-quantification conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/openms/quantification/openms.yaml
# Build: docker build -t openms-quantification -f Dockerfile .
# Save:  docker save openms-quantification -o /home/luosg/Database/env/openms-quantification.tar

FROM continuumio/miniconda3:latest

COPY openms.yaml /tmp/openms.yaml

RUN conda env create -f /tmp/openms.yaml && \
    conda clean -afy && \
    rm /tmp/openms.yaml

ENV CONDA_DEFAULT_ENV=openms-quantification
ENV PATH="/opt/conda/envs/openms-quantification/bin:$PATH"

CMD ["bash"]
