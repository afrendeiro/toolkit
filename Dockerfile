# FROM gitpod/workspace-full
FROM gitpod/workspace-full-vnc

USER root

# Install bedtools
RUN apt-get update \
    && sudo apt-get install -y --no-install-recommends \
        bedtools \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
    && rm -rf /var/lib/apt/lists/*

# Install R and bioconductor libraries
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        r-base \
        r-bioc-deseq2 \
        r-bioc-preprocesscore \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
    && rm -rf /var/lib/apt/lists/*

USER gitpod

ENV PYTHONPATH=/home/gitpod/.local/lib/python3.7/site-packages/

# Install IPython
RUN pip3 install --user ipython

# Install Python dependencies of ngs-toolkit
RUN pip3 install --user -r \
        https://raw.githubusercontent.com/afrendeiro/toolkit/master/requirements/requirements.txt \
    && pip3 install --user -r \
        https://raw.githubusercontent.com/afrendeiro/toolkit/master/requirements/requirements.test.txt \
    && pip3 install --user git+https://github.com/afrendeiro/combat.git

# # Install library
# RUN pip3 install --user \
#     git+https://github.com/afrendeiro/toolkit.git#egg=ngs-toolkit[testing]

USER root

ENV PATH="/home/gitpod/.local/bin:${PATH}"
