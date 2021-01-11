#! /usr/bin/env python

import logging

def get_logger(mod_name):

    logger = logging.getLogger(mod_name)
    logger.propagate = False

    if not logger.hasHandlers(): 
        ch = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m/%d/%Y %H:%M:%S')
        ch.setFormatter(formatter)
        logger.addHandler(ch)

    logger.setLevel(logging.INFO)

    return logger

