import logging
import os
import sys
from logging.handlers import RotatingFileHandler
from logging import FileHandler

def get_log_path():
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "run.log")

def get_logger(name):
    file = get_log_path()
    log_formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(filename)s(line %(lineno)d) %(message)s')

    #file_handler = RotatingFileHandler(file, mode='w+', maxBytes=5*1024*1024,
    #                                 backupCount=5, encoding=None, delay=0)
    file_handler = FileHandler(file, mode='a', encoding=None)
    file_handler.setFormatter(log_formatter)
    file_handler.setLevel(logging.DEBUG)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(log_formatter)
    console_handler.setLevel(logging.INFO)

    #error_handler = logging.StreamHandler()
    #error_handler.setFormatter(log_formatter)
    #error_handler.setLevel(logging.ERROR)

    app_logger = logging.getLogger(name)
    app_logger.setLevel(logging.DEBUG)
    app_logger.addHandler(file_handler)
    app_logger.addHandler(console_handler) #todo zapnout
   # app_logger.addHandler(error_handler)

    return app_logger