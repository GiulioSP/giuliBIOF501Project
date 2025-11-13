FROM continuumio/miniconda3

WORKDIR /app

COPY amrlearn.yml .

RUN conda install -n base -c conda-forge mamba && \
    mamba env create -f amrlearn.yml --retry-clean-cache && \
    conda clean -afy

ARG ENV_NAME=amrlearn

SHELL ["conda", "run", "-n", "amrlearn", "/bin/bash", "-c"]

CMD ["python"]
