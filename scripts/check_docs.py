#!/usr/bin/env python3
"""Documentation coverage check.

Every function definition in the scanned directories must be preceded by a
block comment, per section 3a of the implementation plan. The rule exists
because this codebase is full of quantities that are dimensionally
interchangeable and physically are not -- Hz against rad/s, slant against
ground range, azimuth line rate against transmit PRF, electromagnetic against
acoustic wavelength -- and a comment stating units and conventions is the
cheapest defence available against confusing them.

The check is deliberately crude: it looks for a line that opens a function body
at column zero and walks backwards past any attributes to find a `*/`. Crude is
fine. The point is that the rule is enforced from the first commit rather than
becoming a cleanup task before publication.
"""

import re
import sys
from pathlib import Path

# A function definition: a line starting at column zero that ends in `{` and
# contains a parameter list. Excludes control keywords, struct/enum/union
# definitions, and anything that is plainly a declaration ending in `;`.
DEF_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_ *]*\**\s*[A-Za-z_][A-Za-z0-9_]*\s*\(")
KEYWORDS = {"if", "for", "while", "switch", "return", "else", "do",
            "struct", "union", "enum", "typedef"}


def find_undocumented(path: Path) -> list[tuple[int, str]]:
    """Return (line number, text) for each undocumented function in a file."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    problems = []

    for i, line in enumerate(lines):
        if not DEF_RE.match(line):
            continue
        first_word = line.split("(")[0].split()[0] if line.split("(")[0].split() else ""
        if first_word in KEYWORDS:
            continue
        # A definition opens a body: `{` on this line, or on the next
        # non-blank one for multi-line parameter lists.
        opens_body = line.rstrip().endswith("{")
        if not opens_body:
            for j in range(i + 1, min(i + 4, len(lines))):
                nxt = lines[j].strip()
                if not nxt:
                    continue
                if nxt.endswith("{"):
                    opens_body = True
                if nxt.endswith(";"):
                    break
                break
        if not opens_body:
            continue

        # Walk backwards past blank lines and preprocessor guards to find the
        # end of a block comment.
        k = i - 1
        while k >= 0 and (not lines[k].strip() or lines[k].lstrip().startswith("#")):
            k -= 1
        if k < 0 or not lines[k].rstrip().endswith("*/"):
            problems.append((i + 1, line.strip()[:76]))

    return problems


def main(argv: list[str]) -> int:
    """Scan every .c file under the given directories."""
    if len(argv) < 2:
        print("usage: check_docs.py DIR [DIR ...]", file=sys.stderr)
        return 2

    total = 0
    scanned = 0
    for root in argv[1:]:
        for path in sorted(Path(root).rglob("*.c")):
            scanned += 1
            for lineno, text in find_undocumented(path):
                print(f"{path}:{lineno}: undocumented function: {text}")
                total += 1

    if total:
        print(f"\n{total} function(s) lack the required block comment "
              f"(see IMPLEMENTATION_PLAN.md section 3a)")
        return 1

    print(f"docs coverage: {scanned} file(s) scanned, every function documented")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
