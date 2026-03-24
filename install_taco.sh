#!/usr/bin/env bash
set -euo pipefail

printf '[info] TACO installer\n'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

read -r -p "Install folder [~/opt]: " INSTALL_BASE
INSTALL_BASE="${INSTALL_BASE:-~/opt}"
INSTALL_BASE="${INSTALL_BASE/#\~/$HOME}"

APP_DIR="${INSTALL_BASE}/TACO"
BIN_DIR="${INSTALL_BASE}/bin"
mkdir -p "$APP_DIR" "$BIN_DIR"

find_latest_taco_script() {
  local dir="$1"
  local latest=""

  latest="$(
    find "$dir" -maxdepth 1 -type f -name 'TACO-*.sh' -printf '%f\n' \
      | sed -E 's/^TACO-([0-9]+(\.[0-9]+)*)\.sh$/\1\t&/' \
      | sort -V \
      | tail -n 1 \
      | cut -f2
  )"

  if [[ -z "$latest" && -f "$dir/TACO.sh" ]]; then
    latest="TACO.sh"
  fi

  [[ -n "$latest" ]] && printf '%s\n' "$dir/$latest"
}

PIPELINE_SRC="$(find_latest_taco_script "$SCRIPT_DIR")"
if [[ -z "${PIPELINE_SRC:-}" || ! -f "$PIPELINE_SRC" ]]; then
  echo "[error] No pipeline script found in $SCRIPT_DIR" >&2
  echo "[error] Expected TACO.sh or TACO-<version>.sh" >&2
  exit 1
fi

echo "[info] Using pipeline: $(basename "$PIPELINE_SRC")"

PIPELINE_DST="$APP_DIR/TACO.sh"
cp -f "$PIPELINE_SRC" "$PIPELINE_DST"
chmod +x "$PIPELINE_DST"

ENV_SRC="$SCRIPT_DIR/taco-env.yml"
ENV_DST="$APP_DIR/taco-env.yml"

if [[ ! -f "$ENV_SRC" ]]; then
  echo "[error] taco-env.yml not found next to installer" >&2
  exit 1
fi

cp -f "$ENV_SRC" "$ENV_DST"

INSTALLER_DST="$APP_DIR/install_taco.sh"
SELF_SRC="$SCRIPT_DIR/$(basename "$0")"
cp -f "$SELF_SRC" "$INSTALLER_DST"
chmod +x "$INSTALLER_DST"

ENV_NAME="taco"

# detect solver
if command -v micromamba >/dev/null 2>&1; then
  SOLVER="micromamba"
elif command -v conda >/dev/null 2>&1; then
  SOLVER="conda"
elif command -v mamba >/dev/null 2>&1; then
  SOLVER="mamba"
else
  echo "[error] Need micromamba, conda, or mamba" >&2
  exit 1
fi

echo "[info] Using solver: $SOLVER"

create_env() {
  if [[ "$SOLVER" == "micromamba" ]]; then
    micromamba env create -y -n "$ENV_NAME" -f "$ENV_DST" || \
    micromamba env update -y -n "$ENV_NAME" -f "$ENV_DST"
  elif [[ "$SOLVER" == "conda" ]]; then
    conda env create -y -n "$ENV_NAME" -f "$ENV_DST" || \
    conda env update -y -n "$ENV_NAME" -f "$ENV_DST"
  else
    mamba env create -y -n "$ENV_NAME" -f "$ENV_DST" || \
    mamba env update -y -n "$ENV_NAME" -f "$ENV_DST"
  fi
}

create_env

# install redundans
if ! command -v git >/dev/null 2>&1; then
  echo "[error] git is required to install redundans" >&2
  exit 1
fi

REDUNDANS_DIR="$APP_DIR/redundans"
echo "[info] Installing redundans into $REDUNDANS_DIR"

if [[ ! -d "$REDUNDANS_DIR/.git" ]]; then
  git clone --recursive https://github.com/Gabaldonlab/redundans.git "$REDUNDANS_DIR"
else
  git -C "$REDUNDANS_DIR" pull --ff-only
  git -C "$REDUNDANS_DIR" submodule update --init --recursive
fi

cd "$REDUNDANS_DIR"

if [[ -x "bin/.compile.sh" ]]; then
  bash bin/.compile.sh
elif [[ -f "INSTALL.sh" ]]; then
  bash INSTALL.sh
else
  echo "[error] redundans installer not found" >&2
  exit 1
fi

# wrapper
cat > "$BIN_DIR/redundans.py" <<EOF
#!/usr/bin/env bash
exec "$REDUNDANS_DIR/redundans.py" "\$@"
EOF
chmod +x "$BIN_DIR/redundans.py"

# TACO launcher
cat > "$BIN_DIR/TACO.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
if command -v micromamba >/dev/null 2>&1; then
  exec micromamba run -n $ENV_NAME "$APP_DIR/TACO.sh" "\$@"
elif command -v conda >/dev/null 2>&1; then
  exec conda run -n $ENV_NAME "$APP_DIR/TACO.sh" "\$@"
elif command -v mamba >/dev/null 2>&1; then
  exec mamba run -n $ENV_NAME "$APP_DIR/TACO.sh" "\$@"
else
  echo "[error] No environment runner found" >&2
  exit 127
fi
EOF
chmod +x "$BIN_DIR/TACO.sh"

# PATH
BASHRC="$HOME/.bashrc"
PATH_LINE='export PATH="$HOME/opt/bin:$PATH"'

if ! grep -Fq "$PATH_LINE" "$BASHRC" 2>/dev/null; then
  {
    echo
    echo '# TACO path'
    echo "$PATH_LINE"
  } >> "$BASHRC"
fi

echo
echo "[ok] Installation complete"
echo "[ok] Main app: $APP_DIR/TACO.sh"
echo "[ok] Launcher: $BIN_DIR/TACO.sh"