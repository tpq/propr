FROM bioconductor/bioconductor_docker:RELEASE_3_12
LABEL author="Suzanne Jin" \
      github="suzannejin/propr" \
      source-Dockerfile="https://github.com/suzannejin/propr/blob/master/Dockerfile" \
      description="A Docker container with all the packages required to run R package <propr>"


# install last version propr & dependencies
# this is my dev version
# at some point the official version should be updated too
RUN R -e "BiocManager::install('ALDEx2')"
RUN R -e "install.packages('cccrm')"
RUN R -e "install.packages('compositions')"
RUN R -e "install.packages('corpcor')"
RUN R -e "install.packages('data.table')"
RUN R -e "install.packages('ggdendro')"
RUN R -e "BiocManager::install('limma')"
RUN R -e "install.packages('plotly')"
RUN R -e "install.packages('ppcor')"
RUN R -e "install.packages('reshape2')"
RUN R -e "install.packages('rgl')"
RUN R -e "install.packages('rmarkdown')"
RUN R -e "install.packages('vegan')"
RUN R -e "devtools::install_github('suzannejin/propr@6f1a7fdd477cd4e6d61fee1ab582ff8dbaa1e60b')"


# install graflex & dependencies
RUN R -e "BiocManager::install('AnnotationDbi')"
RUN R -e "devtools::install_github('tpq/graflex')"