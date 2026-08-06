# Auto-generated Dockerfile for UmiTools conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/UmiTools/UmiTools.yaml
# Build: docker build -t UmiTools -f Dockerfile .
# Save:  docker save UmiTools -o /home/luosg/Database/env/UmiTools.tar

FROM continuumio/miniconda3:latest

COPY UmiTools.yaml /tmp/UmiTools.yaml

RUN conda env create -f /tmp/UmiTools.yaml && \
    conda clean -afy && \
    rm /tmp/UmiTools.yaml

ENV CONDA_DEFAULT_ENV=UmiTools
ENV PATH="/opt/conda/envs/UmiTools/bin:$PATH"

CMD ["bash"]
