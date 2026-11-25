# Auto-generated Dockerfile for DR_report conda environment
# Source: workflow/Omics/modules/DR_report/DR_report.yaml
# Build: docker build -t DR_report -f Dockerfile .
# Save:  docker save DR_report -o /data/16T/Database/env/DR_report.tar

FROM continuumio/miniconda3:latest

COPY DR_report.yaml /tmp/DR_report.yaml

RUN conda env create -f /tmp/DR_report.yaml && \
    conda clean -afy && \
    rm /tmp/DR_report.yaml

ENV CONDA_DEFAULT_ENV=DR_report
ENV PATH="/opt/conda/envs/DR_report/bin:$PATH"

CMD ["bash"]
