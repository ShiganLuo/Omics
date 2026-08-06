# Auto-generated Dockerfile for mimseq conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/mimseq/mimseq.yaml
# Build: docker build -t mimseq -f Dockerfile .
# Save:  docker save mimseq -o /home/luosg/Database/env/mimseq.tar

FROM continuumio/miniconda3:latest

COPY mimseq.yaml /tmp/mimseq.yaml

RUN conda env create -f /tmp/mimseq.yaml && \
    conda clean -afy && \
    rm /tmp/mimseq.yaml

ENV CONDA_DEFAULT_ENV=mimseq
ENV PATH="/opt/conda/envs/mimseq/bin:$PATH"

CMD ["bash"]
