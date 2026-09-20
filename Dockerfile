# Container image pinning the full software environment used for the SCART
# spatial transcriptomics analysis.
#
# Build (from the repository root):
#   docker build -t ghcr.io/<ORG>/scart-spatial:v1.0.0 .
#
# On Apple Silicon, force the architecture others will run:
#   docker build --platform linux/amd64 -t ghcr.io/<ORG>/scart-spatial:v1.0.0 .
#
# Run:
#   docker run --rm -v /path/to/data:/data -v $(pwd)/results:/results \
#     -e SCART_DATA=/data ghcr.io/<ORG>/scart-spatial:v1.0.0 \
#     python scripts/regenerate_fig4C.py

# Base image tag: confirm this tag still exists before the first build
#   docker pull condaforge/miniforge3:24.9.2-0
# If it does not, pick a current tag from https://hub.docker.com/r/condaforge/miniforge3/tags
# and pin it here.  Do NOT use :latest - it defeats the point of the image.
FROM condaforge/miniforge3:24.9.2-0

LABEL org.opencontainers.image.title="scart-spatial"
LABEL org.opencontainers.image.description="Analysis environment for the SCART spatial transcriptomics study (murine H22 HCC)"
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.source="https://github.com/<ORG>/SCART-spatial-2026"

SHELL ["/bin/bash", "-lc"]

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg

# Fonts matter: the figures specify Arial and fall back to DejaVu Sans.  Without
# a font package the container silently renders a different face than the
# manuscript figures.
RUN apt-get update && apt-get install -y --no-install-recommends \
        fonts-dejavu-core \
        git \
        procps \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# If conda downloads stall behind a slow link, uncomment the mirror below.
# (Tsinghua TUNA mirror - commonly needed for builds inside mainland China.)
# ---------------------------------------------------------------------------
# RUN conda config --set show_channel_urls yes && \
#     conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge && \
#     conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda && \
#     pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple

WORKDIR /opt/scart

# Copy the environment spec alone first so that edits to the analysis scripts do
# not invalidate the (slow) dependency layer.
COPY environment.yml /opt/scart/environment.yml

RUN mamba env create -f /opt/scart/environment.yml && \
    mamba clean --all --yes && \
    find /opt/conda -follow -type f -name '*.a'    -delete && \
    find /opt/conda -follow -type f -name '*.pyc'  -delete

# Put the environment on PATH so `python` resolves without `conda activate`.
ENV PATH=/opt/conda/envs/scart-spatial/bin:$PATH \
    CONDA_DEFAULT_ENV=scart-spatial

COPY . /opt/scart

# Data is mounted at runtime, never baked into the image (GEO is the source).
ENV SCART_DATA=/data
VOLUME ["/data", "/results"]

# Fail the build rather than ship a broken image.
RUN python -c "import scanpy, anndata, squidpy, decoupler, matplotlib, scipy; \
print('scanpy', scanpy.__version__); \
print('anndata', anndata.__version__); \
print('squidpy', squidpy.__version__); \
print('decoupler', decoupler.__version__)"

CMD ["python", "-c", "import scanpy; scanpy.logging.print_header()"]
