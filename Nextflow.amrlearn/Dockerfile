FROM continuumio/miniconda3

WORKDIR /app

COPY amrlearn.yml /tmp/amrlearn.yml

#RUN conda install -n base -c conda-forge mamba && \
#    mamba env create -f amrlearn.yml --retry-clean-cache && \
#    conda clean -afy

RUN conda env create -f /tmp/amrlearn.yml && \
    conda clean --all -y

ARG ENV_NAME=amrlearn

SHELL ["conda", "run", "-n", "amrlearn", "/bin/bash", "-c"]

CMD ["python"]

WORKDIR /workspace

ENTRYPOINT ["conda", "run", "-n", "amrlearn"]