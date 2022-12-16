import logging

logging.getLogger("matplotlib").setLevel(level=logging.CRITICAL)


def setup_logger():
    logger = logging.getLogger("GASSY")
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler()
    format = logging.Formatter("GASSY:%(levelname)s:%(message)s")
    handler.setFormatter(format)
    logger.addHandler(handler)
    return logger


logger = setup_logger()
