import os
import tempfile

import pytest

from ..utils.param_utils import delete_param, set_param


@pytest.fixture
def example_file():
    text = b"""
%% param_Noh_3d.txt
% parameter file for 3d Noh problem

InitCondFile                          ./IC
ICFormat                              3

OutputDir                             ./output/
SnapshotFileBase                      snap
SnapFormat                            3
NumFilesPerSnapshot                   1
NumFilesWrittenInParallel             1

    """
    with tempfile.NamedTemporaryFile(delete=False) as fp:
        fp.write(text)
        fp.close()
        yield fp.name
        os.remove(fp.name)


def test_replacing_param(example_file):
    old_line = "ICFormat                              3\n"
    param = "ICFormat"
    value = 1
    expected_new_line = "ICFormat                              1\n"

    with open(example_file, "r") as f:
        old_lines = f.readlines()
    assert old_line in old_lines
    old_idx = old_lines.index(old_line)

    set_param(example_file, param=param, value=value, value_column=39)
    with open(example_file, "r") as f:
        new_lines = f.readlines()

    assert expected_new_line in new_lines
    assert old_idx == new_lines.index(expected_new_line)


def test_space_after_long_param_name(example_file):
    long_param_name = "Lorem_ipsum_dolor_sit_amet_consectetur_adipiscing_elit"
    value = 1
    expected_line = long_param_name + " 1\n"

    set_param(example_file, param=long_param_name, value=value)
    with open(example_file) as f:
        lines = f.readlines()

    assert lines[-1] == expected_line


def test_adding_param_after_other(example_file):
    param = "ICStyle"
    value = 2
    expected_line = "ICStyle                               2\n"
    after_param = "ICFormat"

    set_param(example_file, param=param, value=value, after_param="ICFormat")

    with open(example_file) as f:
        lines = f.readlines()
    assert expected_line in lines
    previous_line = lines[lines.index(expected_line) - 1]
    assert after_param in previous_line


def test_adding_param_before_other(example_file):
    param = "ICStyle"
    value = 2
    expected_line = "ICStyle                               2\n"
    before_param = "ICFormat"

    set_param(example_file, param=param, value=value, before_param="ICFormat")

    with open(example_file) as f:
        lines = f.readlines()
    assert expected_line in lines
    next_line = lines[lines.index(expected_line) + 1]
    assert before_param in next_line


def test_removing_param(example_file):
    param = "ICFormat"
    delete_param(example_file, param=param)
    with open(example_file) as f:
        lines = f.readlines()

    matching_lines = [l for l in lines if param in l]
    assert not matching_lines
