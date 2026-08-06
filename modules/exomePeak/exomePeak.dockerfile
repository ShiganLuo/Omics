# Auto-generated Dockerfile for exomePeak conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/exomePeak/exomePeak.yaml
# Build: docker build -t exomePeak -f Dockerfile .
# Save:  docker save exomePeak -o /home/luosg/Database/env/exomePeak.tar

FROM continuumio/miniconda3:latest

COPY exomePeak.yaml /tmp/exomePeak.yaml

RUN conda env create -f /tmp/exomePeak.yaml && \
    conda clean -afy && \
    rm /tmp/exomePeak.yaml

ENV CONDA_DEFAULT_ENV=exomePeak
ENV PATH="/opt/conda/envs/exomePeak/bin:$PATH"

CMD ["bash"]
