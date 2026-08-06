# Auto-generated Dockerfile for scRNAseq_scanpy conda environment
# Source: /home/luosg/Data/genomeStability/workflow/Omics/modules/scanpy/scRNAseq_scanpy.yaml
# Build: docker build -t scRNAseq_scanpy -f Dockerfile .
# Save:  docker save scRNAseq_scanpy -o /home/luosg/Database/env/scRNAseq_scanpy.tar

FROM continuumio/miniconda3:latest

COPY scRNAseq_scanpy.yaml /tmp/scRNAseq_scanpy.yaml

RUN conda env create -f /tmp/scRNAseq_scanpy.yaml && \
    conda clean -afy && \
    rm /tmp/scRNAseq_scanpy.yaml

ENV CONDA_DEFAULT_ENV=scRNAseq_scanpy
ENV PATH="/opt/conda/envs/scRNAseq_scanpy/bin:$PATH"

CMD ["bash"]
