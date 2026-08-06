# Auto-generated Dockerfile for cnvkit conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/cnvkit/cnvkit.yaml
# Build: docker build -t cnvkit -f Dockerfile .
# Save:  docker save cnvkit -o /home/luosg/Database/env/cnvkit.tar

FROM continuumio/miniconda3:latest

COPY cnvkit.yaml /tmp/cnvkit.yaml

RUN conda env create -f /tmp/cnvkit.yaml && \
    conda clean -afy && \
    rm /tmp/cnvkit.yaml

ENV CONDA_DEFAULT_ENV=cnvkit
ENV PATH="/opt/conda/envs/cnvkit/bin:$PATH"

CMD ["bash"]
