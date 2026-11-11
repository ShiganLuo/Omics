# Auto-generated Dockerfile for report conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/RIPseq_report/RIPseq_report.yaml
# Build: docker build -t report -f Dockerfile .
# Save:  docker save report -o /home/luosg/Database/env/report.tar

FROM continuumio/miniconda3:latest

COPY RIPseq_report.yaml /tmp/RIPseq_report.yaml

RUN conda env create -f /tmp/RIPseq_report.yaml && \
    conda clean -afy && \
    rm /tmp/RIPseq_report.yaml

ENV CONDA_DEFAULT_ENV=report
ENV PATH="/opt/conda/envs/report/bin:$PATH"

CMD ["bash"]
