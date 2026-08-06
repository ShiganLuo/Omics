# Auto-generated Dockerfile for karrseq conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/karrseq/karrseq.yaml
# Build: docker build -t karrseq -f Dockerfile .
# Save:  docker save karrseq -o /home/luosg/Database/env/karrseq.tar

FROM continuumio/miniconda3:latest

COPY karrseq.yaml /tmp/karrseq.yaml

RUN conda env create -f /tmp/karrseq.yaml && \
    conda clean -afy && \
    rm /tmp/karrseq.yaml

ENV CONDA_DEFAULT_ENV=karrseq
ENV PATH="/opt/conda/envs/karrseq/bin:$PATH"

CMD ["bash"]
