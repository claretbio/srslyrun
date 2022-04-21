from unittest.mock import patch

import sys
from contextlib import contextmanager
from io import StringIO

from workflow.cli import main


@contextmanager
def capture_cli(argv):
    """This context manager captures the standard output and error of a
    function, as well as setting the command line arguments of
    sys.argv. This is meant for testing functions that are used as a
    command line interface.
    """
    new_out, new_err = StringIO(), StringIO()
    with patch.object(sys, "argv", argv):
        with patch.object(sys, "stdout", new_out):
            with patch.object(sys, "stderr", new_err):
                yield sys.stdout, sys.stderr


def test_main():
    argv = ["foobar-run"]
    with capture_cli(argv) as (stdout, stderr):
        main()
    assert stderr.getvalue() == ""
    assert stdout.getvalue() == "Hello World!\n"


def test_main_message():
    argv = ["foobar-run", "--message", "Goodbye"]
    with capture_cli(argv) as (stdout, stderr):
        main()
    assert stderr.getvalue() == ""
    assert stdout.getvalue() == "Goodbye\n"
