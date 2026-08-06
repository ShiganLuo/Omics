# Auto-generated Dockerfile for telomere conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/telomere/telomere.yaml
# Build: docker build -t telomere -f Dockerfile .
# Save:  docker save telomere -o /home/luosg/Database/env/telomere.tar

FROM continuumio/miniconda3:latest

COPY telomere.yaml /tmp/telomere.yaml

RUN conda env create -f /tmp/telomere.yaml && \
    conda clean -afy && \
    rm /tmp/telomere.yaml

ENV CONDA_DEFAULT_ENV=telomere
ENV PATH="/opt/conda/envs/telomere/bin:$PATH"

CMD ["bash"]
