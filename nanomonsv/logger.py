#! /usr/bin/env python3

import logging
from os import path

global_log_dir: str = ''

def set_logger_dir(log_dir: str) -> None:
    global global_log_dir
    global_log_dir = log_dir

def get_logger(file_name='nano_log.tsv') -> logging.Logger:

    logger: logging.Logger = logging.getLogger("nano_logger")

    if not logger.hasHandlers():

        file_path: str = path.join(global_log_dir, file_name)
        ch = logging.FileHandler(file_path, mode='w')
        formatter = logging.Formatter(
            '%(asctime)s\t%(filename)s\t%(funcName)s\t%(lineno)d\t'
            '%(levelname)s\t%(message)s', datefmt='%m/%d/%Y %H:%M:%S')
        ch.setFormatter(formatter)
        logger.addHandler(ch)

        logger.propagate = False
        logger.setLevel(logging.INFO)

    return logger
