# Auto-generated Dockerfile for seqprep conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/seqprep/seqprep.yaml
# Build: docker build -t seqprep -f Dockerfile .
# Save:  docker save seqprep -o /home/luosg/Database/env/seqprep.tar

FROM continuumio/miniconda3:latest

COPY seqprep.yaml /tmp/seqprep.yaml

RUN conda env create -f /tmp/seqprep.yaml && \
    conda clean -afy && \
    rm /tmp/seqprep.yaml

ENV CONDA_DEFAULT_ENV=seqprep
ENV PATH="/opt/conda/envs/seqprep/bin:$PATH"

CMD ["bash"]
