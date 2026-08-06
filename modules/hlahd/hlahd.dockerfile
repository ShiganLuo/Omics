# Auto-generated Dockerfile for hlahd conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/hlahd/hlahd.yaml
# Build: docker build -t hlahd -f Dockerfile .
# Save:  docker save hlahd -o /home/luosg/Database/env/hlahd.tar

FROM continuumio/miniconda3:latest

COPY hlahd.yaml /tmp/hlahd.yaml

RUN conda env create -f /tmp/hlahd.yaml && \
    conda clean -afy && \
    rm /tmp/hlahd.yaml

ENV CONDA_DEFAULT_ENV=hlahd
ENV PATH="/opt/conda/envs/hlahd/bin:$PATH"

CMD ["bash"]
