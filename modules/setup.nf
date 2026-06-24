/*
 * SETUP_DEPS — one-time machine setup for disorder predictors + bigBedToBed.
 *
 * Runs bin/setup_external_programs.sh which:
 *   • clones AIUPred / AIUPred-Binding from GitHub
 *   • creates conda envs discanvis_aiupred + discanvis_deepcoil
 *   • installs bigBedToBed static binary into $CONDA_PREFIX/bin
 *   • writes .aiupred_python.txt / .deepcoil_python.txt to ext_programs dir
 *
 * storeDir: outputs land in ${ext_programs}/.nf_setup_done so this process
 * is skipped on every subsequent run on the same machine (files already exist).
 *
 * Downstream consumers (DISORDER_MAP, COILEDCOILS_MAP, POLYMORPHISM_MAP)
 * receive the sentinel file as an input to enforce ordering, plus the detected
 * python-path files so they can auto-configure aiupred_python / deepcoil_python
 * without hardcoded machine paths.
 */

process SETUP_DEPS {
    tag   { "setup_deps" }
    label 'process_low'

    storeDir {
        def base = params.ext_programs ?: "${workflow.projectDir}/External_Programs"
        workflow.stubRun
            ? "${base}/.nf_setup_done_stub"
            : "${base}/.nf_setup_done"
    }

    output:
    path 'setup.done',          emit: done
    path 'aiupred_python.txt',  emit: aiupred_python
    path 'deepcoil_python.txt', emit: deepcoil_python

    script:
    def target    = params.ext_programs ?: "${workflow.projectDir}/External_Programs"
    def iupred_env = params.iupred3_url  ? "IUPRED3_URL=${params.iupred3_url}" : ""
    """
    set -euo pipefail

    # Run the setup script (idempotent: skips anything already installed)
    ${iupred_env} bash "${workflow.projectDir}/bin/setup_external_programs.sh" "${target}"

    # Read the python paths written by the setup script; fall back to conda detection
    _aiupred_py=""
    _deepcoil_py=""
    if [[ -f "${target}/.aiupred_python.txt" ]]; then
        _aiupred_py="\$(cat "${target}/.aiupred_python.txt" | tr -d '[:space:]')"
    fi
    if [[ -f "${target}/.deepcoil_python.txt" ]]; then
        _deepcoil_py="\$(cat "${target}/.deepcoil_python.txt" | tr -d '[:space:]')"
    fi

    # If still empty, try to detect from the conda base directly
    if [[ -z "\${_aiupred_py}" ]] && command -v conda >/dev/null 2>&1; then
        _base="\$(conda info --base 2>/dev/null)"
        [[ -x "\${_base}/envs/discanvis_aiupred/bin/python" ]]  && _aiupred_py="\${_base}/envs/discanvis_aiupred/bin/python"
        [[ -x "\${_base}/envs/discanvis_deepcoil/bin/python" ]] && _deepcoil_py="\${_base}/envs/discanvis_deepcoil/bin/python"
    fi

    echo "\${_aiupred_py}"  > aiupred_python.txt
    echo "\${_deepcoil_py}" > deepcoil_python.txt
    touch setup.done

    echo "[SETUP_DEPS] aiupred_python  = \${_aiupred_py:-<not found>}"
    echo "[SETUP_DEPS] deepcoil_python = \${_deepcoil_py:-<not found>}"
    """

    stub:
    """
    echo "" > aiupred_python.txt
    echo "" > deepcoil_python.txt
    touch setup.done
    """
}
