import logging
import os
import sys
from logging import FileHandler

def get_log_path():
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "run.log")

def get_logger(name, file_handler_level=logging.DEBUG, console_handler_level=logging.INFO):
    file = get_log_path()
    log_formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(filename)s(line %(lineno)d) %(message)s')

    file_handler = FileHandler(file, mode='a', encoding=None)
    file_handler.setFormatter(log_formatter)
    file_handler.setLevel(file_handler_level)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(log_formatter)
    console_handler.setLevel(console_handler_level)

    app_logger = logging.getLogger(name)
    app_logger.setLevel(logging.DEBUG)
    app_logger.addHandler(file_handler)
    app_logger.addHandler(console_handler)

    return app_logger