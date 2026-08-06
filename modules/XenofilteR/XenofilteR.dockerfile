# Auto-generated Dockerfile for XenofilteR conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/XenofilteR/XenofilteR.yaml
# Build: docker build -t XenofilteR -f Dockerfile .
# Save:  docker save XenofilteR -o /home/luosg/Database/env/XenofilteR.tar

FROM continuumio/miniconda3:latest

COPY XenofilteR.yaml /tmp/XenofilteR.yaml

RUN conda env create -f /tmp/XenofilteR.yaml && \
    conda clean -afy && \
    rm /tmp/XenofilteR.yaml

ENV CONDA_DEFAULT_ENV=XenofilteR
ENV PATH="/opt/conda/envs/XenofilteR/bin:$PATH"

CMD ["bash"]
