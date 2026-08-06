# Auto-generated Dockerfile for centromere conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/centromere/centromere.yaml
# Build: docker build -t centromere -f Dockerfile .
# Save:  docker save centromere -o /home/luosg/Database/env/centromere.tar

FROM continuumio/miniconda3:latest

COPY centromere.yaml /tmp/centromere.yaml

RUN conda env create -f /tmp/centromere.yaml && \
    conda clean -afy && \
    rm /tmp/centromere.yaml

ENV CONDA_DEFAULT_ENV=centromere
ENV PATH="/opt/conda/envs/centromere/bin:$PATH"

CMD ["bash"]
