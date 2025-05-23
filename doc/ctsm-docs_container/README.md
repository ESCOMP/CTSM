# The ctsm-docs container
This directory and its Dockerfile are used to build a Docker container for building the CTSM documentation. Unless you're a developer working on the container, you probably don't need to care about anything in here.

## Introduction

This Readme tells you how to update the ctsm-docs Docker container if a need to do so arises—for example, adding a Python module that brings new functionality in the build. After you've followed all these instructions, you will probably want to push an update to [doc-builder](https://github.com/ESMCI/doc-builder) that updates `DEFAULT_DOCKER_IMAGE` in [build_commands.py](https://github.com/ESMCI/doc-builder/blob/master/doc_builder/build_commands.py) to point to the new tag.

## Building

If you actually want to build the container, make sure Docker is running. In the Docker Desktop settings, make sure you've enabled the [`continerd` image store](https://docs.docker.com/desktop/features/containerd/), which allows multi-platform builds. Then do:
```shell
docker buildx build --platform linux/amd64,linux/arm64 -t ghcr.io/escomp/ctsm/ctsm-docs .
```

To use your new version for local testing, you'll need to tell doc-builder to use that image. Call `docker images`, which should return something like this:
```shell
REPOSITORY                      TAG          IMAGE ID       CREATED         SIZE
ghcr.io/escomp/ctsm/ctsm-docs   latest       ab51446519a4   3 seconds ago   233MB
...
```

To test, you can tell `build_docs` to use your new version by adding `--docker-image IMAGE_ID` to your call, where in the example above `IMAGE_ID` is `ab51446519a4`.

## Publishing

### Pushing to GitHub Container Registry
If you want to publish the container, you first need a [GitHub Personal Access Token (Classic)](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens#personal-access-tokens-classic) with the `write:packages` permissions. You can see your existing PAT(C)s [here](https://github.com/settings/tokens). If you don't have one with the right permissions, [this link](https://github.com/settings/tokens/new?scopes=write:packages) should start the setup process for you.

Once you have a PAT(C), you can authenticate in your shell session like so:

```shell
   echo YOUR_PERSONAL_ACCESS_TOKEN_CLASSIC | docker login ghcr.io -u YOUR_USERNAME --password-stdin
```
The leading spaces are intended to prevent this command, which contains your secret PAT(C), from being written to your shell's history file. That at least works in bash... sometimes. To be extra safe, in bash you can do `history -c` and it will clear the session's history entirely.

### Tagging
You'll next need to tag the image. Lots of Docker instructions tell you to use the `latest` tag, and Docker may actually do that for you. However, `latest` can lead to support headaches as users think they have the right version but actually don't. Instead, you'll make a new version number incremented from the [previous one](https://github.com/ESCOMP/CTSM/pkgs/container/ctsm%2Fctsm-docs/versions), in the `vX.Y.Z` format.

Copy the relevant image ID (see `docker images` instructions above) and tag it with your version number like so:
```shell
docker tag ab51446519a4 ghcr.io/escomp/ctsm/ctsm-docs:vX.Y.Z
```

Push to the repo:
```shell
docker push ghcr.io/escomp/ctsm/ctsm-docs:vX.Y.Z
```

Then browse to the [container's GitHub page](https://github.com/ESCOMP/CTSM/pkgs/container/ctsm%2Fctsm-docs) to make sure this all worked and the image is public.

### Updating doc-builder
Since you've updated the container, you will probably want to tell [doc-builder](https://github.com/ESMCI/doc-builder) to use the new one. Open a PR where you change the tag (the part after the colon) in the definition of `DEFAULT_DOCKER_IMAGE` in `doc_builder/build_commands.py`. Remember, **use the version number**, not "latest".

## See also

- [GitHub: Working with the container registry](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry)
- [GitHub: Connecting a repository to a package](https://docs.github.com/en/packages/learn-github-packages/connecting-a-repository-to-a-package)