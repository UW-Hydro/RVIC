import pytest
from time import sleep
import multiprocessing
from rvic.core.multi_proc import LoggingPool as Pool
from rvic.core.multi_proc import LogExceptions
from rvic.core.pycompat import pyrange, bytes_type
from rvic.core.variables import Point
import pickle
from rvic.parameters import store_result, results

multiprocessing.log_to_stderr()


@pytest.fixture()
def pool(scope='function'):
    return Pool()


@pytest.fixture()
def point(scope='function'):
    return Point()


def test_create_log_exceptions():
    LogExceptions(callable)


def test_pool_raises(pool):
    def go_exception():
        print(1)
        sleep(0.5)
        raise Exception()
        print(2)
        return 3

    def go_value_error():
        print(1)
        sleep(0.5)
        raise ValueError()
        print(2)
        return 3

    with pytest.raises(Exception):
        for i in pyrange(3):
            pool.apply_async(go_exception)
        pool.close()
        pool.join()

    with pytest.raises(ValueError):
        for i in pyrange(3):
            pool.apply_async(go_value_error)
        pool.close()
        pool.join()


def test_pool_callback(pool):

    class Foo(object):
        '''simple object to hold a few attributes'''

        def __init__(self, i):
            self.cell_id = i
            self.cube = i ** 3

        def __repr__(self):
            return "%s: %s" % (self.cell_id, self.cube)

        def __str__(self):
            return "%s: %s" % (self.cell_id, self.cube)

    def go_foo(i):
        '''create a foo object'''
        sleep(0.5)
        return Foo(i)

    n = 5

    for i in pyrange(n):
        pool.apply_async(go_foo, store_result, i)
    pool.close()
    pool.join()

    assert type(results) == list
    assert len(results) == n
    for foo in results:
        assert type(foo) == Foo
        assert foo.cube == foo.cell_id ** 3


def test_pickle_point(point):
    pickle_string = pickle.dumps(point)
    assert type(pickle_string) == bytes_type
