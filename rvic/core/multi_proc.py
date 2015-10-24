# -*- coding: utf-8 -*-
'''
multi_proc.py
'''

import multiprocessing


def error(*args):
    """Error function"""
    return multiprocessing.get_logger().error(*args)
