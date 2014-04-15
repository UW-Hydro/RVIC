#!/usr/bin/env python
"""RVIC command line testing interface"""

from __future__ import print_function
import argparse
import pytest
import cProfile
import pstats
import StringIO
from rvic import convert, convolution, parameters
from rvic.core.config import read_config


# -------------------------------------------------------------------- #
def main():
    """
    Run RVIC tests
    """
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
        pytest.main('-x unit')
    if any(i in ['all', 'examples'] for i in args.test_set):
        run_examples(args.examples)
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
    num_tests = len(config_dict.keys())

    for i, (test, test_dict) in enumerate(config_dict.iteritems()):
        print("".center(100, '-'))
        print("Starting Test #{0} of {1}: {2}".format(i+1, num_tests,
                                                      test).center(100))
        print("".center(100, '-'))
        pr = cProfile.Profile()
        pr.enable()

        if test_dict['function'] == 'convert':
            convert.convert(test_dict['config_file'])
        elif test_dict['function'] == 'convolution':
            convolution.convolution(test_dict['config_file'])
        elif test_dict['function'] == 'parameters':
            parameters.parameters(test_dict['config_file'])
        else:
            raise ValueError('Unknow function variable: '
                             '{0}'.format(test_dict['function']))

        pr.disable()
        s = StringIO.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()

        print("".center(100, '-'))
        print("Done With Test #{0} of {1}: {2}".format(i+1, num_tests,
                                                       test).center(100))
        print(".....Printing Profile Information.....".center(100))
        print("".center(100, '-'))
        print(s.getvalue())
        print("".center(100, '-'))

    return


# -------------------------------------------------------------------- #
if __name__ == "__main__":
    main()
# -------------------------------------------------------------------- #
