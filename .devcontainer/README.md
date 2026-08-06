# Claude Code Sandbox

Sandboxed container for running Claude Code unattended on this project. Filesystem access is limited to this repo checkout; outbound network is locked to a small allowlist.

This setup is based on files from Anthropic's [claude-code repo](https://github.com/anthropics/claude-code) with the goal of using development containers as described [here](https://code.claude.com/docs/en/devcontainer).

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

## First-time setup

1. Open this repo in VS Code with the "Dev Containers" extension installed.
2. Command Palette → "Dev Containers: Reopen in Container". First build takes a few minutes.
3. In the VS Code terminal, run `claude` and authenticate.
4. Install the Superpowers plugin:

    ```
    /plugin marketplace add obra/superpowers-marketplace
    /plugin install superpowers@superpowers-marketplace
    ```
   The `~/.claude` directory is on a named volume, so this only needs to happen once—it survives rebuilds.

5. (If building documentation) Authenticate with GitHub:
    ```bash
    gh auth login
    gh auth setup-git
    ```
    You may need to give it a Personal Access Token (classic), which you can create by following [these instructions](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens#creating-a-personal-access-token-classic).

## Running unattended

Inside the container:

```bash
claude --dangerously-skip-permissions
```

The container is the safety boundary, so this flag is fine here. NEVER run with `--dangerously-skip-permissions` outside the container!

## Modifying the firewall

Edit `init-firewall.sh`, then rebuild the container (Command Palette → "Dev Containers: Rebuild Container"). The `postStartCommand` re-runs the script on every start, so changes take effect immediately on rebuild.
