import os
import re
import tempfile

import pytest

from ..utils.crud_file_line import update_or_insert_line, matching_pattern_not_in_file


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


@pytest.fixture
def empty_file():
    with tempfile.NamedTemporaryFile(delete=False) as fp:
        fp.close()
        yield fp.name
        os.remove(fp.name)


class TestLineinfile:

    @staticmethod
    def test_file_unchanged(example_file):
        old_line = "ICFormat                              3\n"
        regex = "ICFormat"
        with open(example_file, "r") as f:
            old_lines = f.readlines()
        assert old_line in old_lines

        update_or_insert_line(example_file, line=old_line, regex=regex, linesep="\n")
        with open(example_file, "r") as f:
            new_lines = f.readlines()

        assert old_lines == new_lines

    @staticmethod
    def test_replacing_line(example_file):
        old_line = "ICFormat                              3\n"
        with open(example_file, "r") as f:
            assert old_line in f

        new_line = "ICFormat                              2\n"
        regex = "ICFormat"
        update_or_insert_line(example_file, line=new_line, regex=regex, linesep="\n")
        with open(example_file, "r") as f:
            lines = f.readlines()

        assert old_line not in lines and new_line in lines

    @staticmethod
    def test_passing_in_regex(example_file):
        old_line = "ICFormat                              3\n"
        new_line = "ICFormat                              2\n"
        regex = re.compile("ICFormat")

        update_or_insert_line(example_file, line=new_line, regex=regex, linesep="\n")
        with open(example_file, "r") as f:
            lines = f.readlines()

        assert old_line not in lines and new_line in lines

    @staticmethod
    def test_adding_line_to_end_of_file(example_file):
        regex = "Missing text"
        line = "Missing text to add\n"
        with open(example_file, "r") as f:
            old_lines = f.readlines()
        assert line not in old_lines

        update_or_insert_line(example_file, line=line, regex=regex, linesep="\n")
        with open(example_file, "r") as f:
            lines = f.readlines()

        assert lines[:-2] == old_lines[:-1]
        # old_lines[-1] has had a newline added
        assert lines[-2] == old_lines[-1] + "\n"
        assert lines[-1] == line

    @staticmethod
    def test_adding_line_without_linesep(example_file):
        regex = "Missing text"
        line_without_linesep = "Missing text to add"
        line = line_without_linesep + "\n"

        update_or_insert_line(example_file, line=line_without_linesep, regex=regex, linesep="\n")
        with open(example_file, "r") as f:
            lines = f.readlines()

        assert lines[-1] == line

    @staticmethod
    def test_adding_line_to_end_of_empty_file(empty_file):
        regex = "Missing text"
        line = "Missing text to add\n"
        with open(empty_file, "r") as f:
            old_lines = f.readlines()
        assert not old_lines

        update_or_insert_line(empty_file, line=line, regex=regex, linesep="\n")
        with open(empty_file, "r") as f:
            lines = f.readlines()

        assert lines[-1] == line

    @staticmethod
    def test_insert_after_match(example_file):
        regex = "Missing text"
        line = "MissingParam                          1\n"
        insertafter = "ICFormat"

        update_or_insert_line(
            example_file, line=line, regex=regex, insert_after_match=insertafter, linesep="\n"
        )
        with open(example_file, "r") as f:
            lines = f.readlines()

        assert line in lines
        idx = lines.index(line)
        assert insertafter in lines[idx - 1]

    @staticmethod
    def test_insert_before_match(example_file):
        regex = "Missing text"
        line = "MissingParam                          1\n"
        insertbefore = "ICFormat"

        update_or_insert_line(
            example_file,
            line=line,
            regex=regex,
            insert_before_match=insertbefore,
            linesep="\n",
        )
        with open(example_file, "r") as f:
            lines = f.readlines()

        assert line in lines
        idx = lines.index(line)
        assert insertbefore in lines[idx + 1]

    @staticmethod
    def test_missing_insert_after_match(example_file):
        regex = "Missing text"
        line = "MissingParam                          1\n"
        insertafter = "Also missing"

        update_or_insert_line(
            example_file, line=line, regex=regex, insert_after_match=insertafter, linesep="\n"
        )
        with open(example_file, "r") as f:
            lines = f.readlines()

        assert lines[-1] == line

    @staticmethod
    def test_remove_line(example_file):
        regex = "ICFormat"
        with open(example_file) as f:
            old_lines = f.readlines()
        matched_lines = [l for l in old_lines if regex in l]
        assert matched_lines

        matching_pattern_not_in_file(example_file, regex=regex)
        with open(example_file) as f:
            lines = f.readlines()

        assert matched_lines[0] not in lines
