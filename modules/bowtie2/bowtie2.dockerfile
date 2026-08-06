# Auto-generated Dockerfile for bowtie2 conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/bowtie2/bowtie2.yaml
# Build: docker build -t bowtie2 -f Dockerfile .
# Save:  docker save bowtie2 -o /home/luosg/Database/env/bowtie2.tar

FROM continuumio/miniconda3:latest

COPY bowtie2.yaml /tmp/bowtie2.yaml

RUN conda env create -f /tmp/bowtie2.yaml && \
    conda clean -afy && \
    rm /tmp/bowtie2.yaml

ENV CONDA_DEFAULT_ENV=bowtie2
ENV PATH="/opt/conda/envs/bowtie2/bin:$PATH"

CMD ["bash"]
