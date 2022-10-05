FROM continuumio/miniconda3:4.10.3p1

RUN apt update -y
RUN apt install -y build-essential

COPY DB_QZ_signac.yaml /home/

RUN conda install mamba -c conda-forge -y
RUN mamba env create -f /home/DB_QZ_signac.yaml

COPY install_signac_environment_packages.R /home/
RUN conda run -n DB_QZ_signac Rscript /home/install_signac_environment_packages.R

RUN echo "conda activate DB_QZ_signac" >> ~/.bashrc
