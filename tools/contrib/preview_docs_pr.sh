#!/usr/bin/env bash
set -ueo pipefail

branch_url="$1"

owner="$(echo $branch_url | cut -d/ -f4)"
fork="$(echo $branch_url | cut -d/ -f5)"
ssh_url="git@github.com:${owner}/${fork}.git"

branch="$(echo $branch_url | cut -d/ -f7)"

clone_dir="${owner}.${fork}.${branch}"

cmd="git clone -b ${branch} ${ssh_url} ${clone_dir}"
echo $cmd

exit 0
