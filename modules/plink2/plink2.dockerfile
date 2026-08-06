# Auto-generated Dockerfile for plink2_population conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/plink2/plink2_population.yaml
# Build: docker build -t plink2_population -f Dockerfile .
# Save:  docker save plink2_population -o /home/luosg/Database/env/plink2_population.tar

FROM continuumio/miniconda3:latest

COPY plink2_population.yaml /tmp/plink2_population.yaml

RUN conda env create -f /tmp/plink2_population.yaml && \
    conda clean -afy && \
    rm /tmp/plink2_population.yaml

ENV CONDA_DEFAULT_ENV=plink2_population
ENV PATH="/opt/conda/envs/plink2_population/bin:$PATH"

CMD ["bash"]
