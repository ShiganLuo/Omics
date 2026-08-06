# Auto-generated Dockerfile for frip_score conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/frip_score/frip_score.yaml
# Build: docker build -t frip_score -f Dockerfile .
# Save:  docker save frip_score -o /home/luosg/Database/env/frip_score.tar

FROM continuumio/miniconda3:latest

COPY frip_score.yaml /tmp/frip_score.yaml

RUN conda env create -f /tmp/frip_score.yaml && \
    conda clean -afy && \
    rm /tmp/frip_score.yaml

ENV CONDA_DEFAULT_ENV=frip_score
ENV PATH="/opt/conda/envs/frip_score/bin:$PATH"

CMD ["bash"]
