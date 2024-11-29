"""Transform an arepo parameters file into c++ code

TODO: We will probably abandon this in favour of calling arepo's read_parameter
function.
"""
from enum import Enum
from pathlib import Path
import re

# FIXME: CLIs
THIS_DIR = Path(__file__).parent
infile = THIS_DIR / "cDobbsParamsTemplate.param"


class ParamType(Enum):
    EMPTY = -1
    COMMENT = 0
    STRING = 1
    INTEGER = 2
    FLOAT = 3
    MEMBERED_STRING = 11
    MEMBERED_INTEGER = 12
    MEMBERED_FLOAT = 13


all_text = infile.read_text()


def param_type(line: str):
    if not line.strip():
        return ParamType.EMPTY

    if line.strip().startswith("%"):
        return ParamType.COMMENT

    line_parts = [p for p in line.split() if p]

    if len(line_parts) < 2 or line_parts[1].strip().startswith("%"):
        raise ValueError("Invalid parameter line: ", line)

    param_name, param_value = line_parts[:2]

    # FIXME: MEMBERED

    return param_value_type(param_value)


def param_value_type(value_string) -> ParamType:
    if re.fullmatch(r"-?\d+", value_string):
        return ParamType.INTEGER

    try:
        float(value_string)
        return ParamType.FLOAT
    except ValueError:
        pass


def param_name_indexed(param_name: str) -> bool:
    # TODO:
    return True
