"""
logger.py
Authors:
    -Stephan Meighen-Berger
Logger handling
"""
import logging


class Logger():
    """
    class: Logger
    Class to inheret logger properties
    Parameters:
        -None
    Returns:
        -None
    """

    def __init__(self):
        """
        function: __init__
        Initializes
        Parameters:
            -None
        Returns:
            -None
        """

    @property
    def logger(self):
        """
        function: logging
        Function called to implement the logger
        Parameters:
            -None
        Returns:
            -None
        """
        component = "{}.{}".format(type(self).__module__, type(self).__name__)
        return logging.getLogger(component)
