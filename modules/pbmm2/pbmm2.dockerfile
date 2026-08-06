# Auto-generated Dockerfile for pbmm2 conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/pbmm2/pbmm2.yaml
# Build: docker build -t pbmm2 -f Dockerfile .
# Save:  docker save pbmm2 -o /home/luosg/Database/env/pbmm2.tar

FROM continuumio/miniconda3:latest

COPY pbmm2.yaml /tmp/pbmm2.yaml

RUN conda env create -f /tmp/pbmm2.yaml && \
    conda clean -afy && \
    rm /tmp/pbmm2.yaml

ENV CONDA_DEFAULT_ENV=pbmm2
ENV PATH="/opt/conda/envs/pbmm2/bin:$PATH"

CMD ["bash"]
