# Custom image for the FETCH_DATABASES process: an Immcantation base with the
# sourcerer package installed on top. sourcerer replaces the fetch/ref2igblast
# shell scripts, so the process container needs it alongside igblast + makeblastdb.
#
# BASE defaults to the Immcantation airrflow image, which already ships igblast,
# changeo and BLAST+. Override it (e.g. the changeo_igblast_wget wave image) with
# --build-arg BASE=... if a lighter or newer-Python base is wanted.
ARG BASE=docker.io/immcantation/airrflow:5.1.0
FROM ${BASE}

# sourcerer is copied in from the build context (the workflow checks it out into
# ./sourcerer_src) rather than pip-installed from a registry, so no PyPI release
# or git credentials are needed to build the test image. Swap this for
# `pip install sourcerer` once it is published.
COPY sourcerer_src /opt/sourcerer_src
RUN pip install --no-cache-dir /opt/sourcerer_src

# Fail the build rather than ship an image that cannot fetch: sourcerer needs
# makeblastdb, and the process reports the igblast version.
RUN sourcerer --version \
    && makeblastdb -version | head -1 \
    && igblastn -version | head -1
