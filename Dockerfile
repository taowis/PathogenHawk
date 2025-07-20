FROM continuumio/miniconda3
WORKDIR /app
COPY . /app
RUN conda env create -f environment.yml
SHELL ["/bin/bash", "-c"]
RUN echo "conda activate pathogenhawk" >> ~/.bashrc
ENV PATH /opt/conda/envs/pathogenhawk/bin:$PATH
CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--allow-root"]