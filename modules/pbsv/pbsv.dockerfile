# Auto-generated Dockerfile for pbsv conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/pbsv/pbsv.yaml
# Build: docker build -t pbsv -f Dockerfile .
# Save:  docker save pbsv -o /home/luosg/Database/env/pbsv.tar

FROM continuumio/miniconda3:latest

COPY pbsv.yaml /tmp/pbsv.yaml

RUN conda env create -f /tmp/pbsv.yaml && \
    conda clean -afy && \
    rm /tmp/pbsv.yaml

ENV CONDA_DEFAULT_ENV=pbsv
ENV PATH="/opt/conda/envs/pbsv/bin:$PATH"

CMD ["bash"]
