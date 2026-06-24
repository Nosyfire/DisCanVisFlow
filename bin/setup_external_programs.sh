#!/usr/bin/env bash
#
# setup_external_programs.sh
# ---------------------------------------------------------------------------
# Make the DisCanVis disorder / coiled-coil predictors reproducible on a clean
# machine by:
#   1. Cloning the public AIUPred GitHub repo(s) into External_Programs/
#      (aiupred-caid3 = disorder lib, AIUPred = AIUPred-Binding lib).
#   2. Installing the IUPred3 / ANCHOR2 lib. These are OUR group's own tools
#      (no third-party licence barrier) — fetched from $IUPRED3_URL if set,
#      otherwise expected from the Zenodo archive or the in-repo External_Programs.
#   3. Creating the two conda envs from envs/*.yml (if conda is available).
#   4. Printing the exact params.aiupred_python / params.deepcoil_python paths
#      to set in the Nextflow profile.
#
# Usage:
#   bin/setup_external_programs.sh [TARGET_DIR]
#   IUPRED3_URL=https://…/iupred3.tar.gz  bin/setup_external_programs.sh
#
#   TARGET_DIR  Optional. Where the External_Programs subdirs are created.
#               Default: <repo root>/External_Programs
#   IUPRED3_URL Optional env var. A tarball URL for IUPred3 (our group hosts it);
#               if set, the script downloads + extracts it automatically.
# ---------------------------------------------------------------------------
set -euo pipefail

# --- Resolve paths relative to this script (so it works from any cwd) --------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
ENVS_DIR="${REPO_ROOT}/envs"

TARGET_DIR="${1:-${REPO_ROOT}/External_Programs}"

echo "=============================================================="
echo " DisCanVis External_Programs setup"
echo "=============================================================="
echo " Repo root   : ${REPO_ROOT}"
echo " Target dir  : ${TARGET_DIR}"
echo " Envs dir    : ${ENVS_DIR}"
echo

mkdir -p "${TARGET_DIR}"

# --- Repo URLs ---------------------------------------------------------------
# Official public AIUPred repository (verified June 2026):
#   https://github.com/doszilab/AIUPred
# It bundles both the disorder predictor and the AIUPred-Binding predictor.
# The pipeline (per CLAUDE.md) expects two subdirs:
#   External_Programs/aiupred-caid3  -> disorder: init_models('disorder') + predict()
#   External_Programs/AIUPred        -> binding : init_models('binding') + predict_binding()
# Both are populated from the same upstream repo. Candidate URLs (in case the
# canonical one moves or you need the CAID3-specific fork):
#   - https://github.com/doszilab/AIUPred      (official, current)
#   - https://github.com/iupred/aiupred        (alt org candidate)
#   - https://github.com/MTB-nrbd/AIUPred      (alt fork candidate)
AIUPRED_DISORDER_URL="https://github.com/doszilab/AIUPred.git"
AIUPRED_BINDING_URL="https://github.com/doszilab/AIUPred.git"

clone_if_missing() {
    local url="$1" dest="$2" name="$3"
    if [[ -d "${dest}" ]]; then
        echo "[skip]  ${name}: directory already exists -> ${dest}"
    else
        echo "[clone] ${name}: ${url}"
        echo "        -> ${dest}"
        git clone --depth 1 "${url}" "${dest}"
    fi
}

echo "--- AIUPred libraries (public GitHub) ------------------------"
clone_if_missing "${AIUPRED_DISORDER_URL}" "${TARGET_DIR}/aiupred-caid3" "AIUPred disorder (aiupred-caid3)"
clone_if_missing "${AIUPRED_BINDING_URL}"  "${TARGET_DIR}/AIUPred"        "AIUPred-Binding (AIUPred)"
echo

# --- IUPred3 / ANCHOR2 (our group's tool — no third-party licence) -----------
echo "--- IUPred3 / ANCHOR2 (our tool) -----------------------------"
IUPRED3_DIR="${TARGET_DIR}/iupred3"
IUPRED3_LIB="${IUPRED3_DIR}/iupred3_lib.py"
if [[ -f "${IUPRED3_LIB}" ]]; then
    echo "[ok]    Found ${IUPRED3_LIB}"
elif [[ -n "${IUPRED3_URL:-}" ]]; then
    echo "[fetch] IUPred3 from \$IUPRED3_URL: ${IUPRED3_URL}"
    mkdir -p "${IUPRED3_DIR}"
    tmp_tar="$(mktemp --suffix=.tar.gz)"
    curl -fsSL "${IUPRED3_URL}" -o "${tmp_tar}"
    # extract; tolerate either a top-level iupred3/ folder or bare files
    tar -xzf "${tmp_tar}" -C "${TARGET_DIR}" || tar -xzf "${tmp_tar}" -C "${IUPRED3_DIR}"
    rm -f "${tmp_tar}"
    if [[ -f "${IUPRED3_LIB}" ]]; then
        echo "[ok]    Installed ${IUPRED3_LIB}"
    else
        echo "[WARN]  Extracted but ${IUPRED3_LIB} not found — check the archive layout."
    fi
else
    echo "[WARN]  IUPred3 library NOT found at: ${IUPRED3_LIB}"
    echo
    echo "        IUPred3 / ANCHOR2 are developed by the Dosztányi lab (ELTE)."
    echo "        Their academic licence prohibits redistribution, so it cannot"
    echo "        be bundled in the pipeline repo or uploaded to GitHub."
    echo
    echo "        Get it by ONE of these methods:"
    echo
    echo "          1. Register + download from the official site (free, academic):"
    echo "               https://iupred2a.elte.hu/download"
    echo "             Then extract into External_Programs/:"
    echo "               tar -xzf iupred3.tar.gz -C ${TARGET_DIR}/"
    echo
    echo "          2. If your group hosts a private copy, re-run with:"
    echo "               IUPRED3_URL=https://your-internal-server/iupred3.tar.gz \\"
    echo "                 bash bin/setup_external_programs.sh"
    echo "             Or set  params.iupred3_url = '...'  in nextflow.config."
    echo
    echo "        Continuing without IUPred3 — IUPred/ANCHOR tracks will be empty."
fi
echo

# --- Conda environments ------------------------------------------------------
echo "--- Conda environments ---------------------------------------"
AIUPRED_ENV="discanvis_aiupred"
DEEPCOIL_ENV="discanvis_deepcoil"

create_env_if_missing() {
    local name="$1" yml="$2"
    if [[ ! -f "${yml}" ]]; then
        echo "[WARN]  Env YAML not found: ${yml} (skipping ${name})"
        return 0
    fi
    # Check by directory (robust even when conda env list format differs in subprocesses)
    local conda_base
    conda_base="$(conda info --base 2>/dev/null || true)"
    local env_dir="${conda_base}/envs/${name}"
    if [[ -d "${env_dir}" ]] || conda env list 2>/dev/null | awk '{print $1}' | grep -qx "${name}"; then
        echo "[skip]  conda env '${name}' already exists (${env_dir})"
    else
        echo "[create] conda env '${name}' from ${yml}"
        # Non-fatal: a broken dependency (e.g. spacy build failure in deepcoil) must
        # not kill the whole setup when the env is expected to already exist.
        conda env create -f "${yml}" \
            || echo "[WARN]  conda env create '${name}' failed — if the env already exists this is harmless; otherwise that predictor will be skipped."
        # discanvis_deepcoil post-install: deepcoil→allennlp==0.9.0→spacy<2.2 cannot
        # be pip-built on modern setuptools (Cython chain broken). Workaround:
        # install spacy 2.1.8 from conda-forge pre-built binary, THEN pip install.
        if [[ "${name}" == "discanvis_deepcoil" ]]; then
            local py="${conda_base}/envs/${name}/bin/python"
            if [[ -x "${py}" ]] && ! "${py}" -c "import deepcoil" 2>/dev/null; then
                # spacy 2.1.x and jsonnet both need C compilation on pip;
                # conda-forge has pre-built binaries for Python 3.7.
                echo "[fix]   Installing spacy 2.1.8 + jsonnet from conda-forge (avoids C build) …"
                conda install -n "${name}" -c conda-forge 'spacy=2.1.8' 'jsonnet' -y \
                    || echo "[WARN]  conda pre-build install failed — trying pip anyway."
                echo "[fix]   pip install allennlp==0.9.0 deepcoil …"
                "${py}" -m pip install allennlp==0.9.0 deepcoil \
                    || echo "[WARN]  deepcoil install failed — CoiledCoils track will be empty."
                # tensorflow 2.9.0 (pulled by deepcoil) is incompatible with protobuf ≥ 4.x
                "${py}" -m pip install "protobuf==3.20.3" \
                    || echo "[WARN]  protobuf downgrade failed — deepcoil may crash on import."
            fi
        fi
    fi
}

if command -v conda >/dev/null 2>&1; then
    create_env_if_missing "${AIUPRED_ENV}"  "${ENVS_DIR}/aiupred.yml"
    create_env_if_missing "${DEEPCOIL_ENV}" "${ENVS_DIR}/deepcoil.yml"
else
    echo "[WARN]  conda not found on PATH — skipping env creation."
    echo "        Install the envs manually once conda is available:"
    echo "          conda env create -f ${ENVS_DIR}/aiupred.yml"
    echo "          conda env create -f ${ENVS_DIR}/deepcoil.yml"
fi
echo

# --- bigBedToBed (UCSC binary for polymorphism track) ------------------------
echo "--- bigBedToBed (UCSC) ---------------------------------------"
if command -v bigBedToBed >/dev/null 2>&1; then
    echo "[ok]    bigBedToBed already in PATH: $(which bigBedToBed)"
elif [[ -n "${CONDA_PREFIX:-}" ]] && [[ -f "${CONDA_PREFIX}/bin/bigBedToBed" ]]; then
    echo "[ok]    bigBedToBed already in \$CONDA_PREFIX/bin"
else
    # Detect platform and install the static UCSC binary into the active conda env
    # (or /usr/local/bin as fallback if no conda env is active).
    _dest="${CONDA_PREFIX:-/usr/local}/bin/bigBedToBed"
    _arch="$(uname -m)"
    if [[ "${_arch}" == "x86_64" ]]; then
        _url="https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed"
    elif [[ "${_arch}" == "aarch64" || "${_arch}" == "arm64" ]]; then
        _url="https://hgdownload.soe.ucsc.edu/admin/exe/linux.aarch64/bigBedToBed"
    else
        echo "[WARN]  Unknown arch ${_arch} — skipping bigBedToBed auto-install."
        _url=""
    fi
    if [[ -n "${_url}" ]]; then
        echo "[install] bigBedToBed from UCSC → ${_dest}"
        curl -fsSL "${_url}" -o "${_dest}" && chmod +x "${_dest}" \
            && echo "[ok]    Installed ${_dest}" \
            || echo "[WARN]  curl failed — polymorphism track will be skipped."
    fi
fi
echo

# --- Detect + write python paths (consumed by Nextflow SETUP_DEPS process) ---
if command -v conda >/dev/null 2>&1; then
    CONDA_BASE="$(conda info --base 2>/dev/null)"
    AIUPRED_PY="${CONDA_BASE}/envs/${AIUPRED_ENV}/bin/python"
    DEEPCOIL_PY="${CONDA_BASE}/envs/${DEEPCOIL_ENV}/bin/python"
else
    AIUPRED_PY=""
    DEEPCOIL_PY=""
fi

# Validate paths exist; blank string signals "not available" to the pipeline
[[ -x "${AIUPRED_PY}"  ]] || AIUPRED_PY=""
[[ -x "${DEEPCOIL_PY}" ]] || DEEPCOIL_PY=""

# Write machine-local path files so SETUP_DEPS (Nextflow) can read them
echo "${AIUPRED_PY}"  > "${TARGET_DIR}/.aiupred_python.txt"
echo "${DEEPCOIL_PY}" > "${TARGET_DIR}/.deepcoil_python.txt"

# --- Print the params paths to set -------------------------------------------
echo "=============================================================="
echo " Set these Nextflow params (e.g. in your profile / config):"
echo "=============================================================="
echo "  params.aiupred_python  = \"${AIUPRED_PY:-<not found — check envs/aiupred.yml>}\""
echo "  params.deepcoil_python = \"${DEEPCOIL_PY:-<not found — check envs/deepcoil.yml>}\""
echo
echo "Done."
