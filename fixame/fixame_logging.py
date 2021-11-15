import logging
import os, errno, sys


logger = logging.getLogger()
FixAME_LOGGING_FORMAT = "[%(levelname)s] - %(asctime)s %(message)s"
FixAME_LOG_LEVEL = logging.INFO


def create_fixame_logging_folders(output_dir, create_folders=[]):
    try:
        full_path = os.path.join(output_dir)
        logger.info(f"FixAME output folder - {output_dir}")
        os.makedirs(os.path.join(full_path))

        for folder in create_folders:
            os.makedirs(os.path.join(full_path, folder))

    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        else:
            logger.exception("It wasn't possible to create FixAME output folder.")
            sys.exit()

    return full_path


def fixame_logging(logger, log_dir, log_level=None, filename=None):
    if log_level is None:
        log_level = FixAME_LOG_LEVEL

    logger.setLevel(log_level)

    if filename is None:
        filename = os.path.join(log_dir,"FixAME_log", "FixAME.log")

    log_formatter = logging.Formatter(FixAME_LOGGING_FORMAT)

    if not logger.handlers:
        terminal_handler = logging.StreamHandler()
        file_handler = logging.FileHandler(filename)

        terminal_handler.setLevel(log_level)
        file_handler.setLevel(log_level)

        log_formatter = logging.Formatter(FixAME_LOGGING_FORMAT)
        terminal_handler.setFormatter(log_formatter)
        file_handler.setFormatter(log_formatter)

        logger.addHandler(terminal_handler)
        logger.addHandler(file_handler)

    return logger
