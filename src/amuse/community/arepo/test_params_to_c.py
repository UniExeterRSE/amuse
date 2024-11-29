import pytest

from .params_to_c import param_type, ParamType

@pytest.mark.parametrize("empty_line", ["", " ", "\n", " \n", "\t"])
def test_empty_strings(empty_line):
    ptype = param_type(empty_line)
    expected_type = ParamType.EMPTY
    assert ptype == expected_type

@pytest.mark.parametrize("comment_line", ["%", "% a", " %", " % c", "%%"])
def test_comments(comment_line):
    ptype = param_type(comment_line)
    expected_type = ParamType.COMMENT
    assert ptype == expected_type

@pytest.mark.parametrize("invalid_line", ["a", "b %"])
def test_invalid_line(invalid_line):
    with pytest.raises(ValueError):
        param_type(invalid_line)

@pytest.mark.parametrize("int_line", ["a 0", "b 1", "c 000" " d -1"])
def test_integer(int_line):
    ptype = param_type(int_line)
    expected_type = ParamType.INTEGER
    assert ptype == expected_type

@pytest.mark.parametrize("float_line", ["a 0.", "b 1.", "c .1e5", " d inf"])
def test_float(float_line):
    ptype = param_type(float_line)
    expected_type = ParamType.FLOAT
    assert ptype == expected_type

