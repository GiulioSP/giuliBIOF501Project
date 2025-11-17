    FROM continuumio/miniconda3:latest

    WORKDIR /app

    COPY aligner.yml .

    RUN conda install -c conda-forge mamba -y && \
    mamba env create -f aligner.yml && \
    conda clean -afy

    ARG ENV_NAME=aligner

    ENTRYPOINT ["conda", "run", "-n", "aligner", "--no-capture-output"]

    CMD ["python"]