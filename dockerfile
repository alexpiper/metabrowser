FROM ubuntu:18.04

LABEL "AUTHOR" "Aurélien DJIAN aurelien.djian@irstea.fr"
LABEL "VERSION" "24/07/2019 - 01"

# Update
RUN apt update && apt upgrade -y

# Configuration tzdata
RUN export DEBIAN_FRONTEND=noninteractive
RUN apt install tzdata -y
RUN ln -fs /usr/share/zoneinfo/Europe/Paris /etc/localtime
RUN dpkg-reconfigure --frontend noninteractive tzdata

# Installation de R et des dépendances
RUN apt-get install git software-properties-common curl libcurl4-openssl-dev libssl-dev ca-certificates wget -y
RUN gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN gpg -a --export E298A3A825C0D65DFD57CBB651716619E084DAB9 | apt-key add -

RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
RUN apt install r-base -y

# Installation des dépendances
RUN R -e "install.packages('shiny')"
RUN R -e "install.packages('ggplot2')"
RUN R -e "install.packages('reshape')"
RUN R -e "install.packages('reshape2')"
RUN R -e "install.packages('ape')"
RUN R -e "install.packages('gridExtra')"
RUN R -e "install.packages('readxl')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "install.packages('DT')"
RUN R -e "install.packages('dplyr')"
RUN R -e "install.packages('glue')"
RUN R -e "install.packages('plotly')"
RUN R -e "install.packages('remotes')"
RUN R -e "install.packages('tibble')"
RUN R -e "install.packages('microbenchmark')"
RUN R -e "install.packages('shinycustomloader')"
RUN R -e "install.packages('shinydashboard')"
RUN R -e "BiocManager::install('phyloseq')"
RUN R -e "remotes::install_github('mahendra-mariadassou/phyloseq-extended', ref='dev')"
RUN R -e "remotes::install_github('kassambara/factoextra')"

# Applications
RUN mkdir -p /opt/easy16S
RUN cd /opt/easy16S
RUN git clone https://gitlab.irstea.fr/cedric.midoux/easy16S.git  /opt/easy16S/easy16S
COPY run.R /opt/easy16S/run.R

WORKDIR /opt/easy16S
