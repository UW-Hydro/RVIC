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


def go_exception(i):
    print('go_exception', i, 0)
    sleep(0.5)
    raise Exception('go_exception raised an Exception')
    print('go_exception', i, 1)
    return ('go_exception', i, 2)


def go_value_error(i):
    print('go_value_error', i, 3)
    sleep(0.5)
    raise ValueError('go_value_error raised an ValueError')
    print('go_value_error', i, 4)
    return ('go_value_error', i, 5)


def test_create_log_exceptions():
    LogExceptions(callable)


def test_pool_raises(pool):

    with pytest.raises(Exception):
        for i in pyrange(3):
            status = pool.apply_async(go_exception, [i])
            status.get()
        pool.close()
        pool.join()

    with pytest.raises(ValueError):
        for i in pyrange(3):
            pool.apply_async(go_value_error, [i])
        pool.close()
        pool.join()


def test_pool_callback(pool):
    n = 5

    for i in pyrange(n):
        pool.apply_async(go_foo, [i], callback=store_result)
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
