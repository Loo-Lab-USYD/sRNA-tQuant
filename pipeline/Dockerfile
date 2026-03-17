FROM mambaorg/micromamba:1.5-jammy AS base

# Metadata
LABEL maintainer="Labhund"
LABEL description="tRNA isodecoder quantification pipeline"

# Install all bioinformatics tools via conda
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Activate the environment
ARG MAMBA_DOCKERFILE_ACTIVATE_ENVIRONMENT=base

# Copy pipeline scripts
COPY --chown=$MAMBA_USER:$MAMBA_USER bin/ /opt/pipeline/bin/
ENV PATH="/opt/pipeline/bin:${PATH}"

# Make scripts executable
RUN chmod +x /opt/pipeline/bin/*.py

WORKDIR /data
