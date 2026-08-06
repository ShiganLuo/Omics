# Auto-generated Dockerfile for genome conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/genome/genome.yaml
# Build: docker build -t genome -f Dockerfile .
# Save:  docker save genome -o /home/luosg/Database/env/genome.tar

FROM continuumio/miniconda3:latest

COPY genome.yaml /tmp/genome.yaml

RUN conda env create -f /tmp/genome.yaml && \
    conda clean -afy && \
    rm /tmp/genome.yaml

ENV CONDA_DEFAULT_ENV=genome
ENV PATH="/opt/conda/envs/genome/bin:$PATH"

CMD ["bash"]
