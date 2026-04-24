#!/usr/bin/env bash
set -ueo pipefail

# TODO: Add --dir option for where to put the clone
# TODO: What happens if someone calls this for a clone dir that already exists?
#   1. Check that it's using the right branch; fail if not
#   2. git pull
#   3. git fleximod update
#   4. Rebuild without -clean by default
#   5. But add -c/--clean to this script as an option
# TODO: Check that URL matches what we expect: https://github.com/$owner/$fork/tree/$branch

branch_url="$1"

owner="$(echo $branch_url | cut -d/ -f4)"
fork="$(echo $branch_url | cut -d/ -f5)"
ssh_url="git@github.com:${owner}/${fork}.git"

branch="$(echo $branch_url | cut -d/ -f7)"

# TODO: May not be within an existing git repo; would cause problems with git-fleximod
# TODO: Default to a subdir of $SCRATCH, if available
clone_dir="preview_docs_pr.${owner}.${fork}.${branch}"

# Clone
cmd="git clone -b ${branch} -o ${owner} ${ssh_url} ${clone_dir}"
echo $cmd
$cmd

# TODO: Merge PR branch into target branch before building

# If on Casper, make sure podman is available
set +u
if [[ "$NCAR_HOST" == "casper" ]]; then
    module load podman
fi
set -u

# Build docs
cd "${clone_dir}"
bin/git-fleximod update
bin/git-fleximod update doc-builder
cd doc
./build_docs -b _build -d

# TODO: Print full path of html/ dir
# TODO: Print paths of directly-affected HTML files

exit 0
