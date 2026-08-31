# shellcheck shell=bash
# Shared helpers for pythonscad.org website and playground SFTP deploys.
# Sourced by composite actions and workflow steps; not executed directly.

validate_web_deploy_config() {
  local label="${1:?config label required}"
  local missing=""
  local var

  for var in DEPLOY_SSH_KEY DEPLOY_HOST_KEY DEPLOY_HOST DEPLOY_USER DEPLOY_PATH; do
    if [ -z "${!var:-}" ]; then
      missing="${missing} ${var}"
    fi
  done

  if [ -n "$missing" ]; then
    {
      echo "Error: required ${label} deploy configuration is missing:${missing}"
      echo
      echo "Populate the matching Settings → Secrets and variables → Actions entries."
      echo "Optional port variable defaults to 22 when unset."
    } >&2
    return 1
  fi
}

setup_web_deploy_ssh() {
  mkdir -p ~/.ssh
  # `printf` (not `echo`) so backslash escapes / option-looking leading bytes
  # in the key material are written byte-faithfully.
  printf '%s' "$DEPLOY_SSH_KEY" > ~/.ssh/deploy_key
  chmod 600 ~/.ssh/deploy_key

  if [ -z "${DEPLOY_HOST_KEY:-}" ]; then
    echo "Error: deploy host key pin is unset." >&2
    return 1
  fi
  printf '%s\n' "$DEPLOY_HOST_KEY" >> ~/.ssh/known_hosts
  chmod 644 ~/.ssh/known_hosts
}

web_sftp_rsync() {
  local port="${DEPLOY_PORT:-22}"
  local remote="${DEPLOY_USER}@${DEPLOY_HOST}:${DEPLOY_PATH%/}"
  local delete_flag=()
  local exclude_flag=()
  local remote_dest

  if [ "${RSYNC_DELETE:-true}" = "true" ]; then
    delete_flag=(--delete)
  fi
  if [ -n "${RSYNC_EXCLUDE:-}" ]; then
    exclude_flag=(--exclude="${RSYNC_EXCLUDE}")
  fi

  if [ -n "${DEPLOY_REMOTE_SUBPATH:-}" ]; then
    remote_dest="${remote}/${DEPLOY_REMOTE_SUBPATH#/}"
  else
    remote_dest="${remote}/"
  fi

  rsync -avz "${delete_flag[@]}" "${exclude_flag[@]}" \
    -e "ssh -i ~/.ssh/deploy_key -p ${port}" \
    "${DEPLOY_SOURCE}" \
    "${remote_dest}"
}

validate_website_build() {
  local site_dir="${1:?site directory required}"

  if [ ! -d "$site_dir" ]; then
    echo "Error: mkdocs output directory ${site_dir} is missing." >&2
    return 1
  fi
  if [ ! -f "${site_dir%/}/index.html" ]; then
    echo "Error: ${site_dir%/}/index.html is missing." >&2
    return 1
  fi
}

validate_playground_dir() {
  local playground_dir="${1:?playground directory required}"

  for required in index.html pythonscad.wasm pythonscad.js pythonscad.data; do
    if [ ! -f "${playground_dir%/}/${required}" ]; then
      echo "Error: playground artifact ${required} is missing under ${playground_dir}." >&2
      return 1
    fi
  done
}
