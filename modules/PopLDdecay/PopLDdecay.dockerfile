# Auto-generated Dockerfile for PopLDdecay_population conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/PopLDdecay/PopLDdecay_population.yaml
# Build: docker build -t PopLDdecay_population -f Dockerfile .
# Save:  docker save PopLDdecay_population -o /home/luosg/Database/env/PopLDdecay_population.tar

FROM continuumio/miniconda3:latest

COPY PopLDdecay_population.yaml /tmp/PopLDdecay_population.yaml

RUN conda env create -f /tmp/PopLDdecay_population.yaml && \
    conda clean -afy && \
    rm /tmp/PopLDdecay_population.yaml

ENV CONDA_DEFAULT_ENV=PopLDdecay_population
ENV PATH="/opt/conda/envs/PopLDdecay_population/bin:$PATH"

CMD ["bash"]
