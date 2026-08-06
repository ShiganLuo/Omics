# Auto-generated Dockerfile for disambiguate conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/disambiguate/disambiguate.yaml
# Build: docker build -t disambiguate -f Dockerfile .
# Save:  docker save disambiguate -o /home/luosg/Database/env/disambiguate.tar

FROM continuumio/miniconda3:latest

COPY disambiguate.yaml /tmp/disambiguate.yaml

RUN conda env create -f /tmp/disambiguate.yaml && \
    conda clean -afy && \
    rm /tmp/disambiguate.yaml

ENV CONDA_DEFAULT_ENV=disambiguate
ENV PATH="/opt/conda/envs/disambiguate/bin:$PATH"

CMD ["bash"]
