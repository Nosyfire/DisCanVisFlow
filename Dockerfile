# DisCanVis pipeline — Docker image
#
# Contains:
#   Python 3.11 + pandas, biopython, numpy, tqdm
#   BLAST+ 2.13.0  (blastp, tblastn, makeblastdb, bl2seq)
#   BLAT + pslCDnaFilter + twoBitToFa  (UCSC genome tools)
#
# Build:
#   docker build -t discanvis-pipeline:latest .
#
# Test:
#   docker run --rm discanvis-pipeline:latest blastp -version
#   docker run --rm discanvis-pipeline:latest blat
#   docker run --rm discanvis-pipeline:latest python3 -c "import pandas, Bio; print('OK')"

FROM ubuntu:22.04

# ── Prevent interactive apt prompts ──────────────────────────────────────────
ENV DEBIAN_FRONTEND=noninteractive

# ── System packages ───────────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3.11 \
        python3.11-dev \
        python3-pip \
        curl \
        wget \
        ca-certificates \
        gzip \
        pigz \
    && rm -rf /var/lib/apt/lists/*

# Make python3.11 the default python3
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1 \
    && update-alternatives --install /usr/bin/python  python  /usr/bin/python3.11 1

# ── Python packages ───────────────────────────────────────────────────────────
RUN pip3 install --no-cache-dir \
        pandas>=1.5 \
        biopython>=1.81 \
        numpy>=1.24 \
        tqdm>=4.64

# ── BLAST+ 2.13.0 ─────────────────────────────────────────────────────────────
RUN wget -q https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-linux.tar.gz \
    && tar xf ncbi-blast-2.13.0+-x64-linux.tar.gz -C /opt \
    && rm ncbi-blast-2.13.0+-x64-linux.tar.gz

ENV PATH="/opt/ncbi-blast-2.13.0+/bin:${PATH}"

# ── UCSC genome tools ─────────────────────────────────────────────────────────
# blat          : cDNA-to-genome alignment
# pslCDnaFilter : PSL output filter
# twoBitToFa    : extract sequences from 2bit genome files
RUN wget -q https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat \
         -O /usr/local/bin/blat \
    && wget -q https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslCDnaFilter \
         -O /usr/local/bin/pslCDnaFilter \
    && wget -q https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa \
         -O /usr/local/bin/twoBitToFa \
    && chmod +x /usr/local/bin/blat \
                /usr/local/bin/pslCDnaFilter \
                /usr/local/bin/twoBitToFa

# ── Working directory ─────────────────────────────────────────────────────────
WORKDIR /pipeline

# ── Sanity check (printed at build time, not cached) ─────────────────────────
RUN blastp -version \
    && blat 2>&1 | head -1 || true \
    && python3 -c "import pandas, Bio, numpy, tqdm; print('Python packages OK')"
