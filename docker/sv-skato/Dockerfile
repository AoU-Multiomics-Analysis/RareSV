FROM mambaorg/micromamba:2.3.2-debian12

RUN micromamba install --yes --name base --strict-channel-priority \
    -c conda-forge -c bioconda \
    r-base=4.5.2 \
    r-tidyverse=2.0.0 \
    r-data.table=1.17.8 \
    r-optparse=1.7.5 \
    r-pgenlibr=0.7.0 \
    r-skat=2.2.5 \
    bioconductor-iranges=2.44.0 \
    bioconductor-qvalue=2.42.0 \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}"
COPY R/sv_set_tests.R R/sv_pgen_workflow.R /opt/sv-skato/
