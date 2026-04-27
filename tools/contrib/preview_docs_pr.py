#!/usr/bin/env python

import urllib.request
import urllib.error
import json
import re
import os
import sys
import argparse
import subprocess
import time
import zipfile
import stat
import shutil

# Get default location in which to clone the code
SCRATCH = os.getenv("SCRATCH")
DEFAULT_extraction_dir_BASENAME = (
    "preview_docs_pr.{}.{}.pr-{}"  # repo owner, repo name, PR number
)
DEFAULT_extraction_dir = os.path.join(
    SCRATCH if SCRATCH else "",
    DEFAULT_extraction_dir_BASENAME,
)

INDENT = 4 * " "


def parse_pr_url(url):
    """Extract owner, repo, and PR number from a GitHub PR URL."""
    pattern = r"https://github\.com/([^/]+)/([^/]+)/pull/(\d+)"
    match = re.match(pattern, url.strip())
    if not match:
        raise ValueError(f"Invalid GitHub PR URL: {url}")
    owner, repo, pr_number = match.groups()
    return owner, repo, int(pr_number)


def make_api_call(url, token):
    headers = {
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
    }
    if token:
        headers["Authorization"] = f"Bearer {token}"

    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req) as resp:
        return json.loads(resp.read().decode())


def fetch_pr_info(owner, repo, pr_number, token=None):
    """Fetch PR metadata from the GitHub API."""

    # Get overall PR info
    pr_url = f"https://api.github.com/repos/{owner}/{repo}/pulls/{pr_number}"
    pr_info = make_api_call(pr_url, token)

    # Get files touched by PR
    files_url = pr_url + "/files"
    pr_files = make_api_call(files_url, token)

    return pr_info, pr_files


def fetch_pr_info_with_mergeability(owner, repo, pr_number, token=None, retries=5):
    wait_time = 10  # seconds
    for attempt in range(retries):
        pr, files = fetch_pr_info(owner, repo, pr_number, token)
        if pr["state"] != "open" or pr.get("mergeable") is not None:
            return pr, files
        print(
            f"  Mergeability not yet computed, waiting {wait_time} seconds and retrying"
            f"({attempt + 1}/{retries})..."
        )
        time.sleep(wait_time)
    raise RuntimeError(f"Mergeability still unknown after {retries} retries.")


def pick_ref(pr):
    """
    Choose the best ref to download, with explanation.
    """
    state = pr["state"]
    merge_commit_sha = pr.get("merge_commit_sha")
    mergeable = pr.get("mergeable")

    # If the PR can't be merged, we can't preview what it'd look like after merge.
    # Note that "mergeable" only refers to whether there are Git conflicts; it doesn't "know"
    # anything about PR tests that are passing or failing.
    if not mergeable:
        print(
            "PR branch has merge conflicts, so merged docs can't be previewed. Exiting."
        )
        sys.exit(1)

    # If PR was closed without merging, GitHub supposedly nulls out merge_commit_sha, so that
    # situation would require special handling.
    if state == "closed" and not pr.get("merged_at"):
        raise NotImplementedError("PR closed without merge")

    if not merge_commit_sha:
        raise RuntimeError("How is this possible?")

    return merge_commit_sha


def download_zip(extraction_dir, token, pr_info, merge_commit_sha):
    print("Downloading zip from GitHub...")
    owner, repo, pr_number = pr_info
    zip_url = f"https://api.github.com/repos/{owner}/{repo}/zipball/{merge_commit_sha}"
    headers = {
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
    }
    if token:
        headers["Authorization"] = f"Bearer {token}"
    req = urllib.request.Request(zip_url, headers=headers)
    pr_number = pr_info[2]
    zip_path = os.path.join(extraction_dir, f"pr-{pr_number}.zip")
    try:
        with urllib.request.urlopen(req) as resp, open(zip_path, "wb") as f:
            f.write(resp.read())
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"Failed to download zip from {zip_url}: {e}") from e
    return zip_path


def extract_zip(extraction_dir, zip_path, overwrite):

    try:
        with zipfile.ZipFile(zip_path, "r") as zf:
            # GitHub zips have a single top-level folder like "owner-repo-<sha>/"
            top_level = zf.namelist()[0].split("/")[0]
            extracted_path = os.path.join(extraction_dir, top_level)

            # If it already exists, it would be nice to skip re-extracting. But we would want to do
            # a git reset to make sure it's actually at the right commit, and unfortunately the
            # downloaded ZIP file doesn't include the git history. So instead, we delete the
            # existing one first.
            if os.path.exists(extracted_path):
                if not overwrite:
                    raise FileExistsError(
                        "Clone directory exists; add -o/--overwrite to overwrite:"
                        f"'{extracted_path}'"
                    )
                print(f"Deleting existing clone dir: '{extracted_path}'")
                shutil.rmtree(extracted_path)

            # Unzip
            print("Unzipping...")
            zf.extractall(extraction_dir)
        os.remove(zip_path)
    except:
        os.remove(zip_path)
        raise

    print(f"Code downloaded to: '{extracted_path}'")
    return extracted_path


def download_pr_code(pr_url, extraction_dir, overwrite, token=None):
    """Download code as if the pull request had been merged"""
    # Fetch PR information via GitHub API
    pr_info = parse_pr_url(pr_url)
    pr, files = fetch_pr_info_with_mergeability(*pr_info, token)

    # Get the ref to download
    merge_commit_sha = pick_ref(pr)

    if not os.path.exists(extraction_dir):
        os.makedirs(extraction_dir)

    # Download the zip (follows redirects automatically)
    zip_path = download_zip(extraction_dir, token, pr_info, merge_commit_sha)

    # Unzip
    code_dir = extract_zip(extraction_dir, zip_path, overwrite)

    return code_dir, files


def build_docs(code_dir, files):
    print("Getting submodules...")
    os.chdir(code_dir)
    os.chmod(path := "bin/git-fleximod", os.stat(path).st_mode | stat.S_IXUSR)
    # subprocess.run([path, "update"], check=True)
    # subprocess.run([path, "update", "doc-builder"], check=True)

    # Do an empty commit to avoid possible errors on git lfs fetch
    subprocess.check_call(["git", "init"])
    subprocess.check_call(
        [
            "git",
            "commit",
            "--allow-empty",
            "-m",
            "Empty commit to avoid git lfs problems",
        ]
    )

    print("Building docs...")
    os.chdir("doc")
    os.chmod(path := "./build_docs", os.stat(path).st_mode | stat.S_IXUSR)
    output = []
    with subprocess.Popen(
        [path, "-b", "_build", "-d", "-c"],
        cwd=os.path.join(code_dir, "doc"),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,  # Merge stderr into stdout
        text=True,
        bufsize=1,  # Line buffering
        env={**os.environ, "PYTHONUNBUFFERED": "1"},  # Force unbuffered output
    ) as process:
        for line in process.stdout:
            print(line, end="")  # Print to screen
            output.append(line)  # Save for later
    output = "".join(output)
    if output and "The HTML pages are in" in output:
        print_files_msg(code_dir, files)


def print_files_msg(code_dir, files):
    build_dir = os.path.join(code_dir, "doc", "_build")
    html_dir = os.path.join(build_dir, "html")
    print(f"\nThe updated files are in {html_dir}")
    print("Doc source files directly touched (not deleted) by PR:")
    for f_dict in files:
        f = f_dict["filename"]
        # Slashes here are platform-independent, because f is returned from GitHub API call
        # Skip deleted or otherwise nonexistent files
        full_path = os.path.join(code_dir, f)
        if f_dict["status"] == "removed" or not os.path.exists(full_path):
            continue

        # Skip files not in doc source
        if not f.startswith("doc/source/"):
            continue

        # Get string to print
        f_print = "/".join(f.split("/")[2:])  # Remove leading doc/source/
        root, extension = os.path.splitext(f_print)
        basename = f.split("/")[-1]  # pylint: disable=use-maxsplit-arg
        if extension in [".rst", ".md"]:
            # These types get converted to HTML
            f_print = root + ".html"
            assert os.path.exists(os.path.join(html_dir, f_print))
        elif os.path.exists(os.path.join(html_dir, "_images", basename)):
            # Image files get put in _build/html/_images/
            f_print = os.path.join("_images", basename)
        else:
            f_print = "[NOT SURE WHERE THIS IS BUILT TO] " + f_print
        print(INDENT + f_print)
    print(
        "Note that changes to these or other files may indirectly affect other doc files!"
        " For example, if one of these files is a new/changed image, or a new/changed text"
        " file that's `include`d somewhere, or a text file with an updated label that's cross-"
        "referenced elsewhere. Or if a file outside doc/source/ that's `include`d in a doc file"
        " got changed."
    )


def parse_args():
    """Parse arguments"""
    parser = argparse.ArgumentParser(
        description="Given a GitHub PR URL, download the merged version and build the docs."
    )

    parser.add_argument(
        "pr_url",
        help="GitHub pull request URL",
        type=str,
    )

    extraction_dir_help_default = DEFAULT_extraction_dir.format(
        "REPO_OWNER", "REPO", "PR_NUM"
    )
    parser.add_argument(
        "--dir",
        "--extraction-dir",
        dest="extraction_dir",
        help=(
            "Directory to download and extract code versions to."
            f"Default: {extraction_dir_help_default}"
        ),
        type=str,
        default=None,
    )

    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite existing clone dir, if any.",
        action="store_true",
    )

    args = parser.parse_args(sys.argv[1:])

    # Get clone dir, if not provided
    pr_info = parse_pr_url(args.pr_url)
    if not args.extraction_dir:
        args.extraction_dir = DEFAULT_extraction_dir.format(*pr_info)
    elif os.path.abspath(args.extraction_dir) == os.getcwd():
        args.extraction_dir = os.path.join(
            os.getcwd(), DEFAULT_extraction_dir_BASENAME.format(*pr_info)
        )
    print(f"Will clone to '{args.extraction_dir}'")

    # Check that clone dir parent exists
    args.extraction_dir = os.path.abspath(args.extraction_dir)
    clone_parent = os.path.dirname(args.extraction_dir)
    if not os.path.exists(clone_parent):
        raise NotImplementedError(
            f"Clone parent directory does not exist: '{clone_parent}'"
        )

    # Check that clone dir parent is not (in) a git repo
    result = subprocess.run(
        ["git", "-C", clone_parent, "rev-parse", "--is-inside-work-tree"],
        check=False,
        capture_output=True,
    )
    if not result.returncode:
        raise RuntimeError(
            "Clone parent directory is (in) a git repo, which would cause git-fleximod problems: "
            f"'{clone_parent}'"
        )

    return args


def main():
    args = parse_args()
    code_dir, files = download_pr_code(args.pr_url, args.extraction_dir, args.overwrite)
    build_docs(code_dir, files)


if __name__ == "__main__":
    main()
