# Auto-generated Dockerfile for vcftools_population conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/vcftools/vcftools_population.yaml
# Build: docker build -t vcftools_population -f Dockerfile .
# Save:  docker save vcftools_population -o /home/luosg/Database/env/vcftools_population.tar

FROM continuumio/miniconda3:latest

COPY vcftools_population.yaml /tmp/vcftools_population.yaml

RUN conda env create -f /tmp/vcftools_population.yaml && \
    conda clean -afy && \
    rm /tmp/vcftools_population.yaml

ENV CONDA_DEFAULT_ENV=vcftools_population
ENV PATH="/opt/conda/envs/vcftools_population/bin:$PATH"

CMD ["bash"]
