#!/usr/bin/env python3
"""
PostToolUse hook that monitors for PR pushes and fetches claude-review results.

When a `gh pr create` or `git push` command is detected, this hook:
1. Extracts the PR number from the command output
2. Waits for the claude-review GitHub Action to complete
3. Fetches any review comments left by the action
4. Returns a system message with the review findings
"""

import json
import sys
import subprocess
import re
import time
import os

# Configuration
MAX_WAIT_SECONDS = 300  # 5 minutes max wait for CI
POLL_INTERVAL_SECONDS = 10
REVIEW_CHECK_NAME = "claude-review"


def get_current_pr_number():
    """Get PR number for current branch if one exists."""
    try:
        result = subprocess.run(
            ["gh", "pr", "view", "--json", "number", "-q", ".number"],
            capture_output=True,
            text=True,
            timeout=30
        )
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()
    except Exception:
        pass
    return None


def extract_pr_number(command_output):
    """Extract PR number from gh pr create or git push output."""
    # Match GitHub PR URLs like https://github.com/owner/repo/pull/123
    pr_match = re.search(r'github\.com/[^/]+/[^/]+/pull/(\d+)', command_output)
    if pr_match:
        return pr_match.group(1)
    return None


def wait_for_claude_review(pr_number):
    """Wait for the claude-review check to complete and return its status."""
    start_time = time.time()

    while time.time() - start_time < MAX_WAIT_SECONDS:
        try:
            result = subprocess.run(
                ["gh", "pr", "checks", pr_number, "--json", "name,state,conclusion"],
                capture_output=True,
                text=True,
                timeout=30
            )

            if result.returncode == 0:
                checks = json.loads(result.stdout)
                for check in checks:
                    if check.get("name") == REVIEW_CHECK_NAME:
                        state = check.get("state", "").upper()
                        if state == "COMPLETED":
                            return check.get("conclusion", "unknown")
                        elif state in ("FAILURE", "ERROR", "CANCELLED"):
                            return state.lower()
        except Exception:
            pass

        time.sleep(POLL_INTERVAL_SECONDS)

    return "timeout"


def fetch_review_comments(pr_number):
    """Fetch review comments from the PR."""
    comments = []

    try:
        # Get PR review comments (inline code comments)
        result = subprocess.run(
            ["gh", "api", f"repos/:owner/:repo/pulls/{pr_number}/comments",
             "--jq", '.[] | select(.user.login == "github-actions[bot]" or .user.type == "Bot") | {path: .path, line: .line, body: .body}'],
            capture_output=True,
            text=True,
            timeout=30
        )

        if result.returncode == 0 and result.stdout.strip():
            for line in result.stdout.strip().split('\n'):
                if line.strip():
                    try:
                        comments.append(json.loads(line))
                    except json.JSONDecodeError:
                        pass

        # Also get issue comments (general PR comments)
        result = subprocess.run(
            ["gh", "api", f"repos/:owner/:repo/issues/{pr_number}/comments",
             "--jq", '.[] | select(.user.login == "github-actions[bot]" or .user.type == "Bot") | {body: .body}'],
            capture_output=True,
            text=True,
            timeout=30
        )

        if result.returncode == 0 and result.stdout.strip():
            for line in result.stdout.strip().split('\n'):
                if line.strip():
                    try:
                        comment = json.loads(line)
                        # Mark as general comment
                        comment["type"] = "general"
                        comments.append(comment)
                    except json.JSONDecodeError:
                        pass

    except Exception as e:
        return [], str(e)

    return comments, None


def format_review_summary(comments, pr_number):
    """Format the review comments into a readable summary."""
    if not comments:
        return f"claude-review completed for PR #{pr_number} with no comments."

    lines = [f"claude-review completed for PR #{pr_number}. Found {len(comments)} comment(s):"]
    lines.append("")

    for i, comment in enumerate(comments, 1):
        if comment.get("type") == "general":
            # Truncate long general comments
            body = comment.get("body", "")
            if len(body) > 500:
                body = body[:500] + "... (truncated)"
            lines.append(f"**General Comment:**")
            lines.append(body)
        else:
            path = comment.get("path", "unknown")
            line_num = comment.get("line", "?")
            body = comment.get("body", "")
            if len(body) > 200:
                body = body[:200] + "..."
            lines.append(f"**{path}:{line_num}**")
            lines.append(body)
        lines.append("")

    return "\n".join(lines)


def main():
    # Read input from stdin
    try:
        input_data = json.load(sys.stdin)
    except json.JSONDecodeError:
        sys.exit(0)

    # Only process Bash tool calls
    if input_data.get("tool_name") != "Bash":
        sys.exit(0)

    tool_input = input_data.get("tool_input", {})
    command = tool_input.get("command", "")
    tool_response = input_data.get("tool_response", {})
    output = str(tool_response.get("output", ""))

    # Check if this is a PR-related command
    is_pr_create = "gh pr create" in command
    is_git_push = "git push" in command and "origin" in command

    if not (is_pr_create or is_git_push):
        sys.exit(0)

    # Try to get PR number
    pr_number = extract_pr_number(output)

    # If git push, might need to get PR from current branch
    if not pr_number and is_git_push:
        pr_number = get_current_pr_number()

    if not pr_number:
        # No PR detected, exit silently
        sys.exit(0)

    # Wait for claude-review to complete
    conclusion = wait_for_claude_review(pr_number)

    if conclusion == "timeout":
        print(json.dumps({
            "continue": True,
            "systemMessage": f"PR #{pr_number} pushed. claude-review is still running (timed out waiting). Check status with: gh pr checks {pr_number}"
        }))
        sys.exit(0)

    if conclusion not in ("success", "SUCCESS"):
        print(json.dumps({
            "continue": True,
            "systemMessage": f"PR #{pr_number}: claude-review finished with status '{conclusion}'. Check the workflow for details."
        }))
        sys.exit(0)

    # Fetch and format review comments
    comments, error = fetch_review_comments(pr_number)

    if error:
        print(json.dumps({
            "continue": True,
            "systemMessage": f"PR #{pr_number}: claude-review completed but couldn't fetch comments: {error}"
        }))
        sys.exit(0)

    summary = format_review_summary(comments, pr_number)

    print(json.dumps({
        "continue": True,
        "systemMessage": summary
    }))
    sys.exit(0)


if __name__ == "__main__":
    main()
