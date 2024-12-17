import os
import re
from pathlib import Path
from re import Pattern
from typing import Union

THIS_DIR = Path(__file__).parent


def update_or_insert_line(
    path: Path,
    line: str,
    regex: Union[str, Pattern],
    *,
    insert_after_match: Union[str, Pattern, None] = None,
    insert_before_match: Union[str, Pattern, None] = None,
    linesep=os.linesep,
):
    """Ensure a line exists in a file

    Search through the file from `path` for a match to `regex` and replace that
    line with `line`.  If there is no match, the line is appended to the end of
    the file, or alternatively after or before matches with `insertafter` or
    `insertbefore` respectively.
    """

    if isinstance(path, str):
        path = Path(path)

    if not path.exists():
        raise ValueError("File doesn't exist: %s", str(path))
    if not path.is_file():
        raise ValueError("Path doesn't represent a file: %s", str(path))

    if not line.endswith(linesep):
        line += linesep

    with path.open() as f:
        lines = f.readlines()

    if line in lines:
        # Nothing to do
        return

    # Scan lines for pattern matches
    matched_regex: list[int] = []
    matched_insertafter: list[int] = []
    matched_insertbefore: list[int] = []
    for idx, l in enumerate(lines):
        if re.match(regex, l):
            matched_regex.append(idx)
        if insert_after_match is not None and re.match(insert_after_match, l):
            matched_insertafter.append(idx)
        if insert_before_match is not None and re.match(insert_before_match, l):
            matched_insertbefore.append(idx - 1)
    n_lines = len(lines)

    # Decision tree of where to insert line
    if matched_regex:
        target_line_number = matched_regex[-1]
        replace = True
    elif matched_insertafter:
        target_line_number = matched_insertafter[-1]
        replace = False
    elif matched_insertbefore:
        target_line_number = matched_insertbefore[-1]
        replace = False
    else:
        # No match - add line to the end of the file
        target_line_number = n_lines
        replace = False

    if target_line_number == n_lines:
        # Check last line ends with newline character
        if lines and not lines[-1].endswith(linesep):
            lines[-1] += linesep

    if replace:
        lines[target_line_number] = line
    else:
        lines.insert(target_line_number + 1, line)

    # Write original contents and the new line to new file
    with path.open("w") as f_out:
        for l in lines:
            f_out.write(l)

    return


def matching_pattern_not_in_file(path: Path, regex: Union[str, Pattern]):
    """
    Remove any matching lines from a file
    """
    if isinstance(path, str):
        path = Path(path)

    if not path.exists():
        raise ValueError("File doesn't exist: %s", str(path))
    if not path.is_file():
        raise ValueError("Path doesn't represent a file: %s", str(path))

    with path.open() as f:
        lines = f.readlines()

    new_lines = [l for l in lines if not re.match(regex, l)]

    with path.open("w") as f:
        for l in new_lines:
            f.write(l)

    return


if __name__ == "__main__":
    update_or_insert_line(THIS_DIR / "param.txt", line="ICFormat", regex="ICFormat", linesep="\n")
