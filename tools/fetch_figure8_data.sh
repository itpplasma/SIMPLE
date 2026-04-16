#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${GITLAB_ACCESS_TOKEN:-}" ]]; then
    echo "GITLAB_ACCESS_TOKEN is not set; skipping figure-8 data fetch"
    exit 0
fi

repo_url="${GITLAB_DATA_REPO_URL:-https://oauth2:${GITLAB_ACCESS_TOKEN}@gitlab.tugraz.at/plasma/data.git}"
figure8_root="${SIMPLE_FIGURE8_DATA_ROOT:-$HOME/data/QUASR/SIMPLE/figure8}"
quasr_root="$(cd "$(dirname "$figure8_root")/../.." && pwd)"
branch="${GITHUB_HEAD_REF:-${GITHUB_REF_NAME:-main}}"

tmpdir="$(mktemp -d)"
cleanup() {
    rm -rf "$tmpdir"
}
trap cleanup EXIT

echo "Fetching figure-8 QUASR data on branch ${branch}"
git clone "$repo_url" "$tmpdir/data"
cd "$tmpdir/data"

git fetch origin
if git ls-remote --exit-code --heads origin "$branch" >/dev/null 2>&1; then
    git checkout "$branch"
else
    echo "Branch ${branch} not found in data repo; falling back to main"
    git checkout main
fi

if [[ -f .gitattributes ]] && grep -q 'filter=lfs' .gitattributes; then
    git lfs install --local
    git config lfs.fetchinclude "QUASR/GVEC/figure8,QUASR/SIMPLE/figure8"
    git config lfs.fetchexclude ""
    git lfs pull
fi

for subdir in QUASR/GVEC/figure8 QUASR/SIMPLE/figure8; do
    if [[ ! -d "$subdir" ]]; then
        echo "Data directory ${subdir} not found in data repo; skipping"
        continue
    fi
    mkdir -p "$quasr_root/$(dirname "${subdir#QUASR/}")"
    rsync -a "$subdir"/ "$quasr_root/${subdir#QUASR/}"/
done

echo "Figure-8 data synchronized to $quasr_root"
