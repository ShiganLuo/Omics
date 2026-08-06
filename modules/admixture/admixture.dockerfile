# Auto-generated Dockerfile for admixture_population conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/admixture/admixture_population.yaml
# Build: docker build -t admixture_population -f Dockerfile .
# Save:  docker save admixture_population -o /home/luosg/Database/env/admixture_population.tar

FROM continuumio/miniconda3:latest

COPY admixture_population.yaml /tmp/admixture_population.yaml

RUN conda env create -f /tmp/admixture_population.yaml && \
    conda clean -afy && \
    rm /tmp/admixture_population.yaml

ENV CONDA_DEFAULT_ENV=admixture_population
ENV PATH="/opt/conda/envs/admixture_population/bin:$PATH"

CMD ["bash"]
