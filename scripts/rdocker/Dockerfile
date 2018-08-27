FROM docker.io/r-base

RUN R -e "install.packages(c('dplyr','ggplot2','grid','reshape2'), repos = 'http://cran.us.r-project.org',dependencies = TRUE)"

COPY . /usr/local/src/myscripts
WORKDIR /usr/local/src/myscripts

CMD ["Rscript", "plotGene_docker.R"]
