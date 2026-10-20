"""Format the generated first-order array facade with the project formatter."""

from pathlib import Path
import re
import subprocess
import sys


path = Path(sys.argv[1])
source = path.read_text()
result = subprocess.run(
    ["fprettify", "--stdout", "-c", ".fprettify", str(path)],
    capture_output=True,
    text=True,
    check=True,
)
if re.search(r"^(ERROR|CRITICAL):", result.stderr, re.MULTILINE):
    raise RuntimeError(result.stderr)
pattern = re.compile(
    r"^[ \t]*subroutine fo_eval_first\(.*?^[ \t]*end subroutine fo_eval_first",
    re.MULTILINE | re.DOTALL,
)
formatted = pattern.search(result.stdout)
if formatted is None or pattern.search(source) is None:
    raise ValueError("generated first-order facade is missing")
source = pattern.sub(lambda _: formatted.group(), source, count=1)
source = re.sub(
    r"^[ \t]*include 'fo_scalar_first.inc'$",
    "    include 'fo_scalar_first.inc'",
    source,
    flags=re.MULTILINE,
)
path.write_text(source)
