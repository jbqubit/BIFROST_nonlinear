# CLAUDE.md — Claude Code Operational Preferences

## Resource Usage

- Default to **low effort** model setting.
- Prefer a single targeted subagent over parallel multi-agent exploration.
- Do not spawn more than 1 subagent unless the task clearly spans multiple independent
  areas.
- Use `/plan` mode for non-trivial tasks before executing.

## Project Guidance

See `AGENT.md` for project-specific guidance including test taxonomy, key invariants, and
what requires user authorization. `AGENT.md` is the authoritative shared spec for all
agents on this project.
