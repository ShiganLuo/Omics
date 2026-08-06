# Auto-generated Dockerfile for fragment_size conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/fragment_size/fragment_size.yaml
# Build: docker build -t fragment_size -f Dockerfile .
# Save:  docker save fragment_size -o /home/luosg/Database/env/fragment_size.tar

FROM continuumio/miniconda3:latest

COPY fragment_size.yaml /tmp/fragment_size.yaml

RUN conda env create -f /tmp/fragment_size.yaml && \
    conda clean -afy && \
    rm /tmp/fragment_size.yaml

ENV CONDA_DEFAULT_ENV=fragment_size
ENV PATH="/opt/conda/envs/fragment_size/bin:$PATH"

CMD ["bash"]
