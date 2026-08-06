# Auto-generated Dockerfile for bcftools_population conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/bcftools/bcftools_population.yaml
# Build: docker build -t bcftools_population -f Dockerfile .
# Save:  docker save bcftools_population -o /home/luosg/Database/env/bcftools_population.tar

FROM continuumio/miniconda3:latest

COPY bcftools_population.yaml /tmp/bcftools_population.yaml

RUN conda env create -f /tmp/bcftools_population.yaml && \
    conda clean -afy && \
    rm /tmp/bcftools_population.yaml

ENV CONDA_DEFAULT_ENV=bcftools_population
ENV PATH="/opt/conda/envs/bcftools_population/bin:$PATH"

CMD ["bash"]
