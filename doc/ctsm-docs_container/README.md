# The ctsm-docs container
This directory and its Dockerfile are used to build a container for building the CTSM documentation. Unless you're a developer working on the container, you probably don't need to care about anything in here.

## Introduction

This Readme tells you how to update the ctsm-docs container if a need to do so arises—for example, adding a Python module that brings new functionality in the build. After you've followed all these instructions, you will probably want to push an update to [doc-builder](https://github.com/ESMCI/doc-builder) that updates `DEFAULT_IMAGE` in [build_commands.py](https://github.com/ESMCI/doc-builder/blob/master/doc_builder/build_commands.py) to point to the new tag.

## Building

If you actually want to build the container, you will need to have Podman installed. We previously used Docker for this, but had to move away from it due to licensing issues. Note that these issues have to do specifically with Docker Desktop, which is required to get the Docker Engine on some platforms. The Docker Engine itself is open-source, so our ctsm-docs container testing and publishing workflows are fine to continue using it. We are also fine to use the Docker Engine via Podman for local builds, which is what's described here.

Once you have Podman installed, you can do:
```shell
podman build --no-cache -t ctsm-docs .
```

To use your new version for local testing, you'll need to tell doc-builder to use that image. Call `podman images`, which should return something like this:
```shell
REPOSITORY                TAG            IMAGE ID      CREATED             SIZE
localhost/ctsm-docs       latest         6464f26339bc  22 seconds ago      241 MB
...
```

To test, you can tell `build_docs` to use your new version by adding `--container-image IMAGE_ID` to your call, where in the example above `IMAGE_ID` is `6464f26339bc`.

## Publishing automatically

The `docker-image-build-publish.yml` workflow makes it so that new versions of the workflow will be published to the GitHub Container Registry whenever changes to the container setup are merged to CTSM's `master` branch. This will fail (as will a similar, no-publish workflow that happens on PRs) unless you specify exactly one new version number in the Dockerfile. This version number will be used as a tag that can be referenced by, e.g., doc-builder.

Lots of container instructions tell you to use the `latest` tag, and indeed the workflow will add that tag automatically. However, actually _using_ `latest` can lead to support headaches as users think they have the right version but actually don't. Instead, you'll make a new version number incremented from the [previous one](https://github.com/ESCOMP/CTSM/pkgs/container/ctsm%2Fctsm-docs/versions), in the `vX.Y.Z` format.

Here's where you need to specify the version number in the Dockerfile:
```docker
LABEL org.opencontainers.image.version="vX.Y.Z"
```
The string there can technically be anything as long as (a) it starts with a lowercase `v` and (b) it hasn't yet been used on a published version of the container.

You can check the results of the automatic publication on the [container's GitHub page](https://github.com/ESCOMP/CTSM/pkgs/container/ctsm%2Fctsm-docs).

### Updating doc-builder
After the new version of the container is published, you will probably want to tell [doc-builder](https://github.com/ESMCI/doc-builder) to use the new one. Open a PR where you change the tag (the part after the colon) in the definition of `DEFAULT_IMAGE` in `doc_builder/build_commands.py`. Remember, **use the version number**, not "latest".

## Publishing manually (NOT recommended)

It's vastly preferable to let GitHub build and publish the new repo using the `docker-image-build-publish.yml` workflow as described above. However, if you need to publish manually for some reason, here's how.

### Building the multi-architecture version

When publishing our container, we need to make sure it can run on either arm64 or amd64 processor architecture. This requires a special build process:
```shell
podman manifest create ctsm-docs-manifest
podman build --platform linux/amd64,linux/arm64  --manifest ctsm-docs-manifest .
```

### Pushing to GitHub Container Registry
If you want to publish the container, you first need a [GitHub Personal Access Token (Classic)](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens#personal-access-tokens-classic) with the `write:packages` permissions. You can see your existing PAT(C)s [here](https://github.com/settings/tokens). If you don't have one with the right permissions, [this link](https://github.com/settings/tokens/new?scopes=write:packages) should start the setup process for you.

Once you have a PAT(C), you can authenticate in your shell session like so:

```bash
# bash: Make it so that, in this session, commands with a leading space are not saved to terminal history
export HISTCONTROL=ignoreboth

# Include leading spaces so that your secret PAT(C) isn't included in terminal history
   echo YOUR_PERSONAL_ACCESS_TOKEN_CLASSIC | podman login ghcr.io -u YOUR_USERNAME --password-stdin
```

### Tagging
You'll next need to tag the image. Lots of container instructions tell you to use the `latest` tag, and Podman may actually add that for you. However, `latest` can lead to support headaches as users think they have the right version but actually don't. Instead, you'll make a new version number incremented from the [previous one](https://github.com/ESCOMP/CTSM/pkgs/container/ctsm%2Fctsm-docs/versions), in the `vX.Y.Z` format.

Copy the relevant image ID (see `podman images` instructions above) and tag it with your version number like so:
```shell
podman tag 6464f26339bc ghcr.io/escomp/ctsm/ctsm-docs:vX.Y.Z
```

Push to the repo:
```shell
podman manifest push --all ctsm-docs-manifest ghcr.io/escomp/ctsm/ctsm-docs:vX.Y.Z
```

Then browse to the [container's GitHub page](https://github.com/ESCOMP/CTSM/pkgs/container/ctsm%2Fctsm-docs) to make sure this all worked and the image is public.

### Updating doc-builder
See "Updating doc-builder" in the "Publishing automatically" section above.

## See also

- [GitHub: Working with the container registry](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry)
- [GitHub: Connecting a repository to a package](https://docs.github.com/en/packages/learn-github-packages/connecting-a-repository-to-a-package)