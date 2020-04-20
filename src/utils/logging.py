import logging,sys

def getLogger(name,filename,stdout=False):
    logger = logging.getLogger('laphasing')
    formatter = logging.Formatter('%(asctime)s | %(levelname)s | %(message)s')
    fhandler  = logging.FileHandler(filename=filename,mode='w')
    fhandler.setLevel(logging.DEBUG)
    fhandler.setFormatter(formatter)
    logger.addHandler(fhandler)
    if stdout:
        shandler = logging.StreamHandler(sys.stdout)
        shandler.setLevel(logging.INFO)
        shandler.setFormatter(formatter)
        logger.addHandler(shandler)
    logger.setLevel(logging.DEBUG)
    return logger
