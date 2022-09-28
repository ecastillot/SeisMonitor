# /**
#  * @author [Emmanuel Castillo]
#  * @email [excastillot@unal.edu.co]
#  * @create date 2021-12-22 09:38:11
#  * @modify date 2021-12-22 09:38:11
#  * @desc [description]
#  */

"""
Utils functions used in the SeisMonitor Module
"""

import logging
import os
import sys

logging.basicConfig(level=logging.INFO,
                   format='%(asctime)s [%(levelname)s] [%(name)s] %(message)s',
                   datefmt='%m-%d %H:%M') 

def validate(func, locals):
    for var, test in func.__annotations__.items():
        try:
            test = test.__args__
            _test_msg = " or ".join(map(str,test))
        except:
            test = test
            _test_msg  = test

        if var == "return":
            continue
        value = locals[var]
        msg = f"Error in {func}: {var} argument must be {_test_msg}"
        assert isinstance(value,test),msg

def printlog(levelname,name,msg):
    """
    Parameters:
    -----------
    levelname: str
        logger levelname
        available: "info","warning,"error"
    name: str
        Subject
    msg: str
        Message that you want to print

    """
    logger = logging.getLogger(name)
    if levelname in ("info","information","INFO","Info","INFORMATION"):
        logger.info(msg)
    elif levelname in ("debug","DEBUG","Debug"):
        logger.debug(msg)
    elif levelname in ("warning","Warning","WARNING"):
        logger.warning(msg)
    elif levelname in ("error","ERROR"):
        logger.error(msg)

def isfile(filepath,overwrite=False):
    """
    Parameters:
    filepath: file path will be saved
    Returns:
    Make the directories needed to save the file.
    If the file is already exist, then ask to the user if want to replace it.
    """


    dirpath = os.path.dirname(filepath)
    if os.path.isdir(dirpath ) == False:
        os.makedirs(dirpath)
    else:
        pass
    if overwrite:
        return True

    if os.path.isfile(filepath) == True:
        while True:
            inp = input(f"{filepath} is already created. Dou you want to replace it? (y or n)")
            if inp.upper() == "Y":
                os.remove(filepath)
                return False
            elif inp.upper() == "N":
                return  True
            else:
                pass
    else:
        return False

