import pytest


def go_and_raise():
    print(1)
    raise Exception()
    print(2)


def test_create_log_exceptions():
    LogExceptions(callable)


def test_create_logging_pool():
    n = 2
    p = LoggingPool(processes=n)
    assert hasattr(p, 'apply_async')
    assert p._processes == n
    assert p.has_logging


def test_logging_pool_raises():
    p = LoggingPool(processes=2)

    with pytest.raises(Exception):
        x = p.apply_async(go_and_raise)
        x.get()
        p.close()
        p.join()
