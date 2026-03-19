#!/usr/bin/env bash
set -euo pipefail

printf '[info] TAMP installer\n'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

read -r -p "Install folder [~/opt]: " INSTALL_BASE
INSTALL_BASE="${INSTALL_BASE:-~/opt}"
INSTALL_BASE="${INSTALL_BASE/#\~/$HOME}"

APP_DIR="${INSTALL_BASE}/TAMP"
BIN_DIR="${INSTALL_BASE}/bin"
mkdir -p "$APP_DIR" "$BIN_DIR"

find_latest_tamp_script() {
  local dir="$1"
  local latest=""

  latest="$(
    find "$dir" -maxdepth 1 -type f -name 'TAMP-*.sh' -printf '%f\n' \
      | sed -E 's/^TAMP-([0-9]+(\.[0-9]+)*)\.sh$/\1\t&/' \
      | sort -V \
      | tail -n 1 \
      | cut -f2
  )"

  if [[ -z "$latest" && -f "$dir/TAMP.sh" ]]; then
    latest="TAMP.sh"
  fi

  [[ -n "$latest" ]] && printf '%s\n' "$dir/$latest"
}

PIPELINE_SRC="$(find_latest_tamp_script "$SCRIPT_DIR")"
if [[ -z "${PIPELINE_SRC:-}" || ! -f "$PIPELINE_SRC" ]]; then
  echo "[error] No pipeline script found in $SCRIPT_DIR" >&2
  echo "[error] Expected TAMP.sh or TAMP-<version>.sh" >&2
  exit 1
fi

echo "[info] Using pipeline: $(basename "$PIPELINE_SRC")"

PIPELINE_DST="$APP_DIR/TAMP.sh"
if [[ ! -e "$PIPELINE_DST" || "$(realpath "$PIPELINE_SRC")" != "$(realpath "$PIPELINE_DST")" ]]; then
  cp "$PIPELINE_SRC" "$PIPELINE_DST"
fi
chmod +x "$PIPELINE_DST"

ENV_SRC="$SCRIPT_DIR/tamp-env.yml"
ENV_DST="$APP_DIR/tamp-env.yml"
if [[ ! -f "$ENV_SRC" ]]; then
  echo "[error] tamp-env.yml not found next to installer" >&2
  exit 1
fi
if [[ ! -e "$ENV_DST" || "$(realpath "$ENV_SRC")" != "$(realpath "$ENV_DST")" ]]; then
  cp "$ENV_SRC" "$ENV_DST"
fi

INSTALLER_DST="$APP_DIR/install_tamp.sh"
SELF_SRC="$SCRIPT_DIR/$(basename "$0")"
if [[ -f "$SELF_SRC" ]]; then
  if [[ ! -e "$INSTALLER_DST" || "$(realpath "$SELF_SRC")" != "$(realpath "$INSTALLER_DST")" ]]; then
    cp "$SELF_SRC" "$INSTALLER_DST"
  fi
  chmod +x "$INSTALLER_DST"
fi

ENV_NAME="tamp"
SOLVER=""
RUNNER=()

if command -v micromamba >/dev/null 2>&1; then
  SOLVER="micromamba"
  RUNNER=(micromamba run -n "$ENV_NAME")
elif command -v conda >/dev/null 2>&1; then
  SOLVER="conda"
  RUNNER=(conda run -n "$ENV_NAME")
elif command -v mamba >/dev/null 2>&1; then
  SOLVER="mamba"
  RUNNER=(mamba run -n "$ENV_NAME")
else
  echo "[error] Need micromamba, conda, or mamba" >&2
  exit 1
fi

echo "[info] Using solver: $SOLVER"

create_env() {
  if [[ "$SOLVER" == "micromamba" ]]; then
    if micromamba env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
      echo "[info] Environment '$ENV_NAME' already exists; updating from $ENV_DST"
      micromamba env update -y -n "$ENV_NAME" -f "$ENV_DST"
    else
      micromamba create -y -n "$ENV_NAME" -f "$ENV_DST"
    fi
  elif [[ "$SOLVER" == "conda" ]]; then
    if conda env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
      echo "[info] Environment '$ENV_NAME' already exists; updating from $ENV_DST"
      conda env update -y -n "$ENV_NAME" -f "$ENV_DST"
    else
      conda env create -y -n "$ENV_NAME" -f "$ENV_DST"
    fi
  else
    if mamba env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
      echo "[info] Environment '$ENV_NAME' already exists; updating from $ENV_DST"
      mamba env update -y -n "$ENV_NAME" -f "$ENV_DST"
    else
      mamba env create -y -n "$ENV_NAME" -f "$ENV_DST"
    fi
  fi
}

create_env

if ! command -v git >/dev/null 2>&1; then
  echo "[error] git is required to install redundans" >&2
  exit 1
fi

REDUNDANS_DIR="$APP_DIR/redundans"
echo "[info] Installing redundans manually into $APP_DIR/redundans"

REDUNDANS_DIR="$APP_DIR/redundans"

if [[ ! -d "$REDUNDANS_DIR/.git" ]]; then
  git clone --recursive https://github.com/Gabaldonlab/redundans.git "$REDUNDANS_DIR"
else
  echo "[info] redundans repo already exists; updating"
  git -C "$REDUNDANS_DIR" pull --ff-only
  git -C "$REDUNDANS_DIR" submodule update --init --recursive
fi

cd "$REDUNDANS_DIR"

# Build bundled third-party tools the way upstream recommends
if [[ -x "bin/.compile.sh" ]]; then
  bash bin/.compile.sh
elif [[ -f "INSTALL.sh" ]]; then
  bash INSTALL.sh
else
  echo "[error] redundans installer files not found (bin/.compile.sh or INSTALL.sh)." >&2
  exit 1
fi

# Ensure launcher exists
mkdir -p "$BIN_DIR"
cat > "$BIN_DIR/redundans.py" <<EOF
#!/usr/bin/env bash
exec "$REDUNDANS_DIR/redundans.py" "\$@"
EOF
chmod +x "$BIN_DIR/redundans.py"

echo "[ok] redundans installed at $REDUNDANS_DIR"

cat > "$BIN_DIR/TAMP.sh" <<EOF2
#!/usr/bin/env bash
set -euo pipefail
if command -v micromamba >/dev/null 2>&1; then
  exec micromamba run -n $ENV_NAME "$APP_DIR/TAMP.sh" "\$@"
elif command -v conda >/dev/null 2>&1; then
  exec conda run -n $ENV_NAME "$APP_DIR/TAMP.sh" "\$@"
elif command -v mamba >/dev/null 2>&1; then
  exec mamba run -n $ENV_NAME "$APP_DIR/TAMP.sh" "\$@"
else
  echo "[error] No environment runner found for TAMP.sh" >&2
  exit 127
fi
EOF2
chmod +x "$BIN_DIR/TAMP.sh"

BASHRC="$HOME/.bashrc"
PATH_LINE='export PATH="$HOME/opt/bin:$PATH"'
if ! grep -Fq "$PATH_LINE" "$BASHRC" 2>/dev/null; then
  {
    echo
    echo '# TAMP path'
    echo "$PATH_LINE"
  } >> "$BASHRC"
  echo "[info] PATH added to ~/.bashrc"
fi

echo
echo "[ok] Installation complete"
echo "[ok] Main app: $APP_DIR/TAMP.sh"
echo "[ok] Launcher: $BIN_DIR/TAMP.sh"
echo "[ok] Redundans wrapper: $BIN_DIR/redundans.py"
echo "[info] Open a new shell or run: source ~/.bashrc"
