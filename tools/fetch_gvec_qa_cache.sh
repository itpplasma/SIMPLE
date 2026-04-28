#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${GITLAB_ACCESS_TOKEN:-}" ]]; then
    echo "GITLAB_ACCESS_TOKEN is not set; skipping GVEC QA cache fetch"
    exit 0
fi

repo_url="${GITLAB_DATA_REPO_URL:-https://oauth2:${GITLAB_ACCESS_TOKEN}@gitlab.tugraz.at/plasma/data.git}"
repo_subdir="${GITLAB_DATA_SIMPLE_GVEC_SUBDIR:-TESTS/SIMPLE/gvec_qa_roundtrip}"
cache_root="${SIMPLE_GVEC_QA_CACHE_ROOT:-$HOME/data/SIMPLE/gvec_qa_roundtrip}"
branch="${GITHUB_HEAD_REF:-${GITHUB_REF_NAME:-main}}"

tmpdir="$(mktemp -d)"
cleanup() {
    rm -rf "$tmpdir"
}
trap cleanup EXIT

echo "Fetching GVEC QA cache from ${repo_subdir} on branch ${branch}"
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
    git config lfs.fetchinclude "$repo_subdir"
    git config lfs.fetchexclude ""
    git lfs pull
    git lfs checkout "$repo_subdir" || true
fi

if [[ ! -d "$repo_subdir" ]]; then
    echo "Cache directory ${repo_subdir} not found in data repo; skipping"
    exit 0
fi

mkdir -p "$cache_root"
rsync -a "$repo_subdir"/ "$cache_root"/
echo "GVEC QA cache synchronized to $cache_root"
