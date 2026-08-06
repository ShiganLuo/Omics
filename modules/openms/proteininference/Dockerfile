# Auto-generated Dockerfile for openms-proteininference conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/openms/proteininference/openms.yaml
# Build: docker build -t openms-proteininference -f Dockerfile .
# Save:  docker save openms-proteininference -o /home/luosg/Database/env/openms-proteininference.tar

FROM continuumio/miniconda3:latest

COPY openms.yaml /tmp/openms.yaml

RUN conda env create -f /tmp/openms.yaml && \
    conda clean -afy && \
    rm /tmp/openms.yaml

ENV CONDA_DEFAULT_ENV=openms-proteininference
ENV PATH="/opt/conda/envs/openms-proteininference/bin:$PATH"

CMD ["bash"]
