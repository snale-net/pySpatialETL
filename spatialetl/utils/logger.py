import logging as log

class myLogger(log.Logger):

    TIMING = 200
    INFO = log.INFO
    DEBUG = log.DEBUG
    WARN = log.WARN
    ERROR = log.ERROR
    CRITICAL = log.CRITICAL

    def __init__(self,name):
        log.Logger.__init__(self,name)
        log.addLevelName(self.TIMING, 'TIMING')
        log.basicConfig(format='[%(levelname)s] %(message)s')

    def debug(self, msg, *args, **kwargs):
        if self.isEnabledFor(self.DEBUG):
            self._log(self.DEBUG, msg, args, **kwargs)

    def info(self, msg, *args, **kwargs):

        if self.isEnabledFor(self.INFO):
            self._log(self.INFO, msg, args, **kwargs)

    def warning(self, msg, *args, **kwargs):
        if self.isEnabledFor(self.WARN):
            self._log(self.WARN, msg, args, **kwargs)

    def timing(self, msg, *args, **kwargs):
        if self.isEnabledFor(self.TIMING):
            self._log(self.TIMING, msg, args, **kwargs)

log.setLoggerClass(myLogger)
logging = log.getLogger("main")








