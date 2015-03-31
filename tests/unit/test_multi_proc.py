from rvic.core.multi_proc import LogExceptions


def test_create_log_exceptions():
    LogExceptions(callable)

# error and logging pool not easily testing in unit test.
