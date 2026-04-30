# CLAUDE.md — Claude Code Operational Preferences

## Resource Usage

- Default to **low effort** model setting.
- Prefer a single targeted subagent over parallel multi-agent exploration.
- Do not spawn more than 1 subagent unless the task clearly spans multiple independent
  areas.
- Use `/plan` mode for non-trivial tasks before executing.

## Dependencies
<!-- The @FILE.md forces claude to read these files into context. -->
@ARCHITECTURE.md
@AGENTS.md
