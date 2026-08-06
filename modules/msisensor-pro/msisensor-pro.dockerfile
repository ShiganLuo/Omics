# Auto-generated Dockerfile for samtools conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/msisensor-pro/msisensro_pro.yaml
# Build: docker build -t samtools -f Dockerfile .
# Save:  docker save samtools -o /home/luosg/Database/env/samtools.tar

FROM continuumio/miniconda3:latest

COPY msisensro_pro.yaml /tmp/msisensro_pro.yaml

RUN conda env create -f /tmp/msisensro_pro.yaml && \
    conda clean -afy && \
    rm /tmp/msisensro_pro.yaml

ENV CONDA_DEFAULT_ENV=samtools
ENV PATH="/opt/conda/envs/samtools/bin:$PATH"

CMD ["bash"]
