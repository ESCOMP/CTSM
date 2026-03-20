#!/usr/bin/env bash
set -e

# Check that documentation links work

# The next line returns code 200 only if the page exists
link_code="$(curl -s -o /dev/null -w "%{http_code}" $URL)"

# Find bad links
set +e
good_links="200"
if [[ "${link_code}" != good_links ]]; then
    echo "One or more links in the documentation failed:" >&2
    echo $link_code  >&2
    exit 1
fi

exit 0
