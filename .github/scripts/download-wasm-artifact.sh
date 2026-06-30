#!/usr/bin/env bash

set -euo pipefail

TARGET_SHA="${1:?target commit sha is required}"
DESTINATION="${2:-wasm-bundle}"
WORKFLOW="${3:-build-wasm.yml}"
ARTIFACT="${4:-pythonscad-wasm-web}"
MAX_WAIT_SECONDS="${WASM_ARTIFACT_MAX_WAIT_SECONDS:-7200}"
POLL_SECONDS="${WASM_ARTIFACT_POLL_SECONDS:-60}"

deadline=$((SECONDS + MAX_WAIT_SECONDS))
start_seconds=${SECONDS}

echo "Waiting for ${WORKFLOW} artifact ${ARTIFACT} from ${TARGET_SHA}"

while true; do
  elapsed_seconds=$((SECONDS - start_seconds))
  run_count="$(gh run list \
    --workflow "${WORKFLOW}" \
    --commit "${TARGET_SHA}" \
    --limit 20 \
    --json databaseId,status,conclusion,createdAt \
    --jq 'length')"

  if [ "${run_count}" = "0" ]; then
    echo "No matching ${WORKFLOW} run exists for ${TARGET_SHA}; cannot download ${ARTIFACT}." >&2
    exit 1
  fi

  run_info="$(gh run list \
    --workflow "${WORKFLOW}" \
    --commit "${TARGET_SHA}" \
    --limit 20 \
    --json databaseId,status,conclusion,createdAt \
    --jq 'if ([.[] | select(.conclusion != "cancelled")] | length) == 0 then empty else ([.[] | select(.conclusion != "cancelled")] | sort_by(.createdAt) | reverse | .[0] | [.databaseId, .status, (.conclusion // "")] | @tsv) end')"

  if [ -n "${run_info}" ]; then
    IFS=$'\t' read -r run_id status conclusion <<<"${run_info}"

    echo "Matched WASM workflow run ${run_id}: status=${status} conclusion=${conclusion:-pending} elapsed=${elapsed_seconds}s"

    if [ "${status}" = "completed" ]; then
      if [ "${conclusion}" != "success" ]; then
        echo "WASM workflow run ${run_id} finished with conclusion=${conclusion}" >&2
        exit 1
      fi
      rm -rf "${DESTINATION}"
      mkdir -p "${DESTINATION}"
      gh run download "${run_id}" --name "${ARTIFACT}" --dir "${DESTINATION}"
      echo "Downloaded ${ARTIFACT} after waiting ${elapsed_seconds}s."
      ls -lh "${DESTINATION}"
      exit 0
    fi
  else
    echo "No non-cancelled matching ${WORKFLOW} run exists for ${TARGET_SHA}; cannot download ${ARTIFACT}." >&2
    exit 1
  fi

  if [ "${SECONDS}" -ge "${deadline}" ]; then
    echo "Timed out waiting for WASM workflow artifact for ${TARGET_SHA}" >&2
    exit 1
  fi

  sleep "${POLL_SECONDS}"
done
