# Auto-generated Dockerfile for report conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/PeakCalling_report/report.yaml
# Build: docker build -t report -f Dockerfile .
# Save:  docker save report -o /home/luosg/Database/env/report.tar

FROM continuumio/miniconda3:latest

COPY report.yaml /tmp/report.yaml

RUN conda env create -f /tmp/report.yaml && \
    conda clean -afy && \
    rm /tmp/report.yaml

ENV CONDA_DEFAULT_ENV=report
ENV PATH="/opt/conda/envs/report/bin:$PATH"

CMD ["bash"]
