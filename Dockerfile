# FROM gitpod/workspace-full-vnc
FROM gitpod/workspace-full

ENV R_BASE_VERSION 3.6.1

USER root

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        apt-utils \
        ed \
        less \
        locales \
        vim-tiny \
        wget \
        ca-certificates \
        fonts-texgyre \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get update

# Use Debian unstable via pinning -- new style via APT::Default-Release
RUN echo "deb http://http.debian.net/debian sid main" > /etc/apt/sources.list.d/debian-unstable.list \
    && echo 'APT::Default-Release "testing";' > /etc/apt/apt.conf.d/default

# Add keys for R packages
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 7638D0442B90D010 \
    && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 04EE7237B7D453EC

# Now install R
RUN apt-get update \
    && apt-get install -t unstable -y --no-install-recommends \
        r-base=${R_BASE_VERSION}-* \
        r-recommended=${R_BASE_VERSION}-* \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
    && rm -rf /var/lib/apt/lists/*

# Install specific bedtools version
RUN wget http://ftp.debian.org/debian/pool/main/b/bedtools/bedtools_2.21.0-1_amd64.deb \
    && dpkg -i bedtools_2.21.0-1_amd64.deb \
    && rm bedtools_2.21.0-1_amd64.deb

RUN R -e "install.packages('BiocManager', repos='https://CRAN.R-project.org')" \
    && R -e "library('BiocManager'); BiocManager::install(c('DESeq2', 'preprocessCore', 'cqn'), version = '3.9')"

RUN pip3 install git+https://github.com/afrendeiro/combat.git \
    && pip3 install ipython

# RUN pip3 install git+https://github.com/afrendeiro/toolkit.git#egg=ngs-toolkit[testing] --user
