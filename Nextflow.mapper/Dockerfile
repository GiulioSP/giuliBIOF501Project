    FROM continuumio/miniconda3:latest

    WORKDIR /app

    COPY mapper.yml .

    RUN conda install -c conda-forge mamba -y && \
    mamba env create -f mapper.yml && \
    conda clean -afy

    ARG ENV_NAME=mapper

    ENTRYPOINT ["conda", "run", "-n", "mapper", "--no-capture-output"]

    CMD ["python"]