from pathlib import Path
from typing import Optional, Union

from .crud_file_line import update_or_insert_line, matching_pattern_not_in_file


def set_param(
    param_file: Path,
    param: str,
    value: Union[int, float, str, bool],
    *,
    after_param: Optional[str] = None,
    before_param: Optional[str] = None,
    value_column: int = 39,
    comment: Optional[str] = None,
    comment_chr: str = "%",
    comment_column: int = 50,
    linesep="\n",
):
    """Set a parameter in a param.txt file, replacing the value if present

    The user supplies the parameter file path, name of the parameter and the
    value.  The function modifies the file, either replacing the line (if the
    parameter is already present) or inserting it otherwise.

    >>> set_param("param.txt", "ICFormat", 2)

    A new parameter line will be added to the end of the file, unless
    `after_param` or `before_param` are supplied (the former taking precedent).

    >>> set_param("param.txt", "MinEgySpec", 0.0, before_param="InitialGasTemp")
    >>> set_param("param.txt", "MinGasTemp", 0.0, after_param="MinEgySpec")

    By default, the format of the modified line will be left justified and
    padded with spaces (at least one) before the value (aligned to
    `value_column`).

    An inline comment can be appended to the parameter line with the `comment`
    argument and, optionally, `comment_chr` and `comment_column`.

    >>> set_param("param.txt", "UnitVelocity_in_cm_per_s", 1e5, value_column=41,
                  comment=" 1 km/sec", comment_chr="%", comment_column=56)
    """
    target_line = "{:{}s} {}".format(param, value_column - 2, value)
    if comment is not None:
        target_line = "{:{}s} {} {}".format(
            target_line, comment_column - 2, comment_chr, comment
        )

    update_or_insert_line(
        param_file,
        regex=param,
        line=target_line,
        insert_after_match=after_param,
        insert_before_match=before_param,
        linesep=linesep,
    )
    return


def delete_param(param_file: Path, param: str):
    """Remove matching line from parameter file"""
    matching_pattern_not_in_file(param_file, regex=param)
    return
