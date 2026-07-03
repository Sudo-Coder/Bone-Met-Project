#!/usr/bin/env bash
# PreToolUse hook for Claude Code (scProj).
# Auto-approves Bash tool calls so the user isn't prompted for every command,
# EXCEPT obviously destructive ones, which are deferred (no decision) so the
# existing deny rules / normal permission prompt still apply.
#
# Wired in via .claude/settings.local.json -> hooks.PreToolUse (matcher "Bash").
set -euo pipefail

input="$(cat)"
tool="$(printf '%s' "$input" | jq -r '.tool_name // empty')"
cmd="$(printf '%s'  "$input" | jq -r '.tool_input.command // empty')"

# Only handle Bash; defer everything else.
[ "$tool" = "Bash" ] || exit 0

lc="$(printf '%s' "$cmd" | tr '[:upper:]' '[:lower:]')"
danger=0

# fork bomb
printf '%s' "$lc" | grep -Eq ':\(\)[[:space:]]*\{[[:space:]]*:\|:&[[:space:]]*\}[[:space:]]*;[[:space:]]*:' && danger=1
# rm with recursive/force flags targeting /, ~, $HOME, or the project root
printf '%s' "$lc" | grep -Eq 'rm[[:space:]].*-[a-z]*[rf][a-z]*[[:space:]].*(/|/\*|~|\$home|/autofs/projects-t3/hussain/scproj)([[:space:]]|/\*|$)' && danger=1
# filesystem / disk destroyers
printf '%s' "$lc" | grep -Eq '(^|[^a-z])(mkfs|shred|fdisk|parted)([^a-z]|$)' && danger=1
printf '%s' "$lc" | grep -Eq 'dd[[:space:]].*of=/dev/' && danger=1
printf '%s' "$lc" | grep -Eq '>[[:space:]]*/dev/(sd|nvme|disk|hd)' && danger=1
# recursive chmod/chown at filesystem root
printf '%s' "$lc" | grep -Eq 'ch(mod|own)[[:space:]]+-[a-z]*r[a-z]*[[:space:]].*[[:space:]]/([[:space:]]|$)' && danger=1
# privilege escalation
printf '%s' "$lc" | grep -Eq '(^|[^a-z])sudo([^a-z]|$)' && danger=1

if [ "$danger" = "1" ]; then
  # Defer: emit nothing so deny rules / the normal prompt decide.
  exit 0
fi

jq -n '{hookSpecificOutput:{hookEventName:"PreToolUse",permissionDecision:"allow",permissionDecisionReason:"Auto-approved Bash in scProj (non-destructive)."}}'
