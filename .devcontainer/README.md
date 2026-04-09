2026-04-09: Copied files here from Anthropic's [claude-code repo](https://github.com/anthropics/claude-code) with the goal of using development containers as described [here](https://code.claude.com/docs/en/devcontainer).

Permalinks to file versions copied:
- [devcontainer.json](https://github.com/anthropics/claude-code/blob/c5600e0b1e9bb6ddf750cf7441c4d4fffbb7c917/.devcontainer/devcontainer.json)
- [Dockerfile](https://github.com/anthropics/claude-code/blob/c5600e0b1e9bb6ddf750cf7441c4d4fffbb7c917/.devcontainer/Dockerfile)
- [init_firewall.sh](https://github.com/anthropics/claude-code/blob/c5600e0b1e9bb6ddf750cf7441c4d4fffbb7c917/.devcontainer/init-firewall.sh)

## CTSM customizations

The vanilla Anthropic setup has been customized for CTSM Python development:

- **Dockerfile**: Installs Miniconda (supports x86_64 and aarch64). Replaces zsh/powerline10k with bash as the default shell.
- **setup-conda-env.sh**: `postCreateCommand` script that creates the `ctsm_pylib` conda environment using `py_env_create` and auto-activates it in new terminal sessions.
- **devcontainer.json**: Configures VS Code with Python extensions (pylint, black-formatter) and settings matching the project's `make lint` and `make black` configurations.
- **init-firewall.sh**: Whitelists conda-forge and PyPI domains so additional packages can be installed after the firewall activates.
