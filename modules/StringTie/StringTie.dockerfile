# Auto-generated Dockerfile for StringTie conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/StringTie/StringTie.yaml
# Build: docker build -t StringTie -f Dockerfile .
# Save:  docker save StringTie -o /home/luosg/Database/env/StringTie.tar

FROM continuumio/miniconda3:latest

COPY StringTie.yaml /tmp/StringTie.yaml

RUN conda env create -f /tmp/StringTie.yaml && \
    conda clean -afy && \
    rm /tmp/StringTie.yaml

ENV CONDA_DEFAULT_ENV=StringTie
ENV PATH="/opt/conda/envs/StringTie/bin:$PATH"

CMD ["bash"]
