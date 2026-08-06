# Auto-generated Dockerfile for openms-psmrescoring conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/openms/psmrescoring/openms.yaml
# Build: docker build -t openms-psmrescoring -f Dockerfile .
# Save:  docker save openms-psmrescoring -o /home/luosg/Database/env/openms-psmrescoring.tar

FROM continuumio/miniconda3:latest

COPY openms.yaml /tmp/openms.yaml

RUN conda env create -f /tmp/openms.yaml && \
    conda clean -afy && \
    rm /tmp/openms.yaml

ENV CONDA_DEFAULT_ENV=openms-psmrescoring
ENV PATH="/opt/conda/envs/openms-psmrescoring/bin:$PATH"

CMD ["bash"]
