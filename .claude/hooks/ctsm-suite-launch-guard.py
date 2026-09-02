#!/usr/bin/env python3
"""PreToolUse guard: ask before a Bash command launches a CTSM system test suite.

A suite may submit hundreds of jobs to a shared queue, billed to the allocation
of whoever launched it, and run for hours on a schedule that person owns. That
decision belongs to a human, so this hook turns a suite launch into a permission
prompt rather than letting it through on an allow rule.

Reads the PreToolUse payload on stdin; prints an "ask" decision when the command
launches a suite, and nothing otherwise (which leaves normal permission handling
in place). See .claude/skills/ctsm-system-tests/SKILL.md.
"""

import json
import os
import shlex
import sys

# Launchers that run a whole suite. create_test is conditional - see below.
SUITE_LAUNCHERS = {"run_sys_tests", "run_sys_tests.py"}
CREATE_TEST = {"create_test", "create_test.py"}

# Words that precede the real command without being it.
WRAPPERS = {
    "qcmd",
    "nohup",
    "time",
    "env",
    "sudo",
    "timeout",
    "xargs",
    "srun",
    "mpirun",
    "mpiexec",
    "exec",
    "command",
    "bash",
    "sh",
    "zsh",
    "python",
    "python3",
    "conda",
    "mamba",
}

SEPARATORS = {"&&", "||", ";", "|", "&", "(", ")", "{", "}", "\n"}

REASON = (
    "This looks like a CTSM system test suite launch ({what}). A suite can submit "
    "hundreds of jobs to the shared queue, billed to the allocation of whoever "
    "launches it, and run for hours - so it is the human's call, not the "
    "assistant's. Assistant-side verification stops at the case build and the "
    "pFUnit unit tests. Approve only if you meant to launch this yourself."
)


def simple_commands(tokens):
    """Split a token list on shell separators into simple commands."""
    current = []
    for token in tokens:
        if token in SEPARATORS:
            if current:
                yield current
            current = []
        else:
            current.append(token)
    if current:
        yield current


def executable_and_args(tokens):
    """Strip wrappers and env assignments; return (basename, remaining args)."""
    index = 0
    while index < len(tokens):
        token = tokens[index]
        if token == "--":
            index += 1
            continue
        if "=" in token and not token.startswith("-") and "/" not in token.split("=")[0]:
            # A VAR=value assignment sitting in front of the command.
            index += 1
            continue
        base = os.path.basename(token)
        if base in WRAPPERS:
            index += 1
            # `python -m ctsm.run_sys_tests` names the target as a module.
            if base.startswith("python") and index < len(tokens) and tokens[index] == "-m":
                index += 1
                if index < len(tokens):
                    return tokens[index].split(".")[-1], tokens[index + 1 :]
                return "", []
            continue
        return base, tokens[index + 1 :]
    return "", []


def what_it_launches(command):
    """Name the suite launch this command performs, or None."""
    try:
        tokens = shlex.split(command, comments=True)
    except ValueError:
        # Unbalanced quotes - fall back to a substring check rather than
        # letting a malformed command slip past the guard.
        if "run_sys_tests" in command:
            return "run_sys_tests"
        if "create_test" in command and "--xml-category" in command:
            return "create_test --xml-category"
        return None

    for tokens in simple_commands(tokens):
        name, args = executable_and_args(tokens)
        if name in SUITE_LAUNCHERS:
            return name
        if name in CREATE_TEST:
            # A single named test is ordinary work; a category is a suite.
            for arg in args:
                if arg == "--xml-category" or arg.startswith("--xml-category="):
                    return "create_test --xml-category"
    return None


def main():
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0

    command = (payload.get("tool_input") or {}).get("command")
    if not isinstance(command, str) or not command.strip():
        return 0

    what = what_it_launches(command)
    if what is None:
        return 0

    print(
        json.dumps(
            {
                "hookSpecificOutput": {
                    "hookEventName": "PreToolUse",
                    "permissionDecision": "ask",
                    "permissionDecisionReason": REASON.format(what=what),
                }
            }
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
