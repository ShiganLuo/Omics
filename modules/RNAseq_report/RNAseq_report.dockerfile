# Auto-generated Dockerfile for RNAseq_report conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/RNAseq_report/RNAseq_report.yaml
# Build: docker build -t RNAseq_report -f Dockerfile .
# Save:  docker save RNAseq_report -o /home/luosg/Database/env/RNAseq_report.tar

FROM continuumio/miniconda3:latest

COPY RNAseq_report.yaml /tmp/RNAseq_report.yaml

RUN conda env create -f /tmp/RNAseq_report.yaml && \
    conda clean -afy && \
    rm /tmp/RNAseq_report.yaml

ENV CONDA_DEFAULT_ENV=RNAseq_report
ENV PATH="/opt/conda/envs/RNAseq_report/bin:$PATH"

CMD ["bash"]
