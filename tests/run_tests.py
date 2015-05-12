#!/usr/bin/env python
"""RVIC command line testing interface"""

import os
import sys
import textwrap
import argparse
import pytest
import cProfile
import pstats
import io
from rvic.core.pycompat import OrderedDict
from rvic import convert, convolution, parameters
from rvic.core.config import read_config
from rvic.core.pycompat import PY3, iteritems

if not os.environ.get('RVIC_TEST_DIR'):
    print('\n$RVIC_TEST_DIR not set.')
    os.environ["RVIC_TEST_DIR"] = os.path.abspath(os.path.dirname(__file__))
    print('Setting to run_tests.py dir: '
          '{0}\n'.format(os.environ["RVIC_TEST_DIR"]))
if not os.environ.get('WORKDIR'):
    print('\n$WORKDIR not set.')
    os.environ["WORKDIR"] = os.environ["RVIC_TEST_DIR"]
    print('Setting to output run_tests.py dir to $WORKDIR: '
          '{0}\n'.format(os.environ["WORKDIR"]))


# -------------------------------------------------------------------- #
def main():
    """
    Run RVIC tests
    """
    exit_code = 0
    # Parse arguments
    parser = argparse.ArgumentParser(description='Test script for RVIC')

    parser.add_argument("test_set", type=str,
                        help="Test set to run",
                        choices=['all', 'unit', 'examples'],
                        default=['all'], nargs='+')
    parser.add_argument("--examples", type=str,
                        help="examples configuration file",
                        default='examples/examples.cfg')
    args = parser.parse_args()

    print('Running Test Set: {0}'.format(args.test_set))

    if any(i in ['all', 'unit'] for i in args.test_set):
        # run unit tests
        exit_code = pytest.main('-x unit')
    if any(i in ['all', 'examples'] for i in args.test_set):
        exit_code = run_examples(args.examples)
    sys.exit(exit_code)
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def run_examples(config_file):
    """ Run examples from config file """
    # ---------------------------------------------------------------- #
    # Read Configuration files
    config_dict = read_config(config_file)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # run tests
    num_tests = len(list(config_dict.keys()))

    test_outcomes = OrderedDict()

    for i, (test, test_dict) in enumerate(iteritems(config_dict)):
        print("".center(100, '-'))
        print("Starting Test #{0} of {1}: {2}".format(i + 1, num_tests,
                                                      test).center(100))
        desc = textwrap.fill(", ".join(test_dict['description']), 100)
        print("Description: {0}".format(desc))
        print("".center(100, '-'))

        if 'processors' in test_dict:
            numofproc = test_dict['processors']
        else:
            numofproc = 1

        pr = cProfile.Profile()
        pr.enable()

        if test_dict['function'] == 'convert':
            try:
                convert.convert(test_dict['config_file'])
                test_outcomes[test] = 'Passed'
            except Exception as e:
                print('Error in conversion example: {0}'.format(test))
                test_outcomes[test] = 'Failed: {0}'.format(e)
        elif test_dict['function'] == 'convolution':
            try:
                convolution.convolution(test_dict['config_file'])
                test_outcomes[test] = 'Passed'
            except Exception as e:
                print('Error in convolution example: {0}'.format(test))
                test_outcomes[test] = 'Failed: {0}'.format(e)
        elif test_dict['function'] == 'parameters':
            try:
                parameters.parameters(test_dict['config_file'],
                                      numofproc=numofproc)
                test_outcomes[test] = 'Passed'
            except Exception as e:
                print('Error in parameters example: {0}'.format(test))
                test_outcomes[test] = 'Failed: {0}'.format(e)
        else:
            raise ValueError('Unknow function variable: '
                             '{0}'.format(test_dict['function']))

        pr.disable()
        if PY3:
            s = io.StringIO()
        else:
            s = io.BytesIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()

        print("".center(100, '-'))
        print("Done With Test #{0} of {1}: {2}".format(i + 1, num_tests,
                                                       test).center(100))
        print(".....Printing Profile Information.....".center(100))
        print("".center(100, '-'))
        print(s.getvalue())
        print("".center(100, '-'))

    print('Done with examples...printing summary')
    for test, outcome in iteritems(test_outcomes):
        print('\t{0}: {1}'.format(test, outcome))

    return 0


# -------------------------------------------------------------------- #
if __name__ == "__main__":
    main()
# -------------------------------------------------------------------- #
