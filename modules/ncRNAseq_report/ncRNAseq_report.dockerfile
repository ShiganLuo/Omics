# Auto-generated Dockerfile for ncRNAseq_report conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/ncRNAseq_report/ncRNAseq_report.yaml
# Build: docker build -t ncRNAseq_report -f Dockerfile .
# Save:  docker save ncRNAseq_report -o /home/luosg/Database/env/ncRNAseq_report.tar

FROM continuumio/miniconda3:latest

COPY ncRNAseq_report.yaml /tmp/ncRNAseq_report.yaml

RUN conda env create -f /tmp/ncRNAseq_report.yaml && \
    conda clean -afy && \
    rm /tmp/ncRNAseq_report.yaml

ENV CONDA_DEFAULT_ENV=ncRNAseq_report
ENV PATH="/opt/conda/envs/ncRNAseq_report/bin:$PATH"

CMD ["bash"]
