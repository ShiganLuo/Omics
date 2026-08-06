# Auto-generated Dockerfile for openms-psmfdr conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/openms/psmfdr/openms.yaml
# Build: docker build -t openms-psmfdr -f Dockerfile .
# Save:  docker save openms-psmfdr -o /home/luosg/Database/env/openms-psmfdr.tar

FROM continuumio/miniconda3:latest

COPY openms.yaml /tmp/openms.yaml

RUN conda env create -f /tmp/openms.yaml && \
    conda clean -afy && \
    rm /tmp/openms.yaml

ENV CONDA_DEFAULT_ENV=openms-psmfdr
ENV PATH="/opt/conda/envs/openms-psmfdr/bin:$PATH"

CMD ["bash"]
