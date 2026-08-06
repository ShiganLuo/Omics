# Auto-generated Dockerfile for featureCounts conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/featureCounts/featureCounts.yaml
# Build: docker build -t featureCounts -f Dockerfile .
# Save:  docker save featureCounts -o /home/luosg/Database/env/featureCounts.tar

FROM continuumio/miniconda3:latest

COPY featureCounts.yaml /tmp/featureCounts.yaml

RUN conda env create -f /tmp/featureCounts.yaml && \
    conda clean -afy && \
    rm /tmp/featureCounts.yaml

ENV CONDA_DEFAULT_ENV=featureCounts
ENV PATH="/opt/conda/envs/featureCounts/bin:$PATH"

CMD ["bash"]
