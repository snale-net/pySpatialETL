import logging as log

class myLogger(log.Logger):

    TIMING = 200
    RUN = 25
    INFO = log.INFO
    DEBUG = log.DEBUG
    WARNING = log.WARN
    ERROR = log.ERROR
    CRITICAL = log.CRITICAL

    def __init__(self,name):
        log.Logger.__init__(self,name)
        log.addLevelName(self.TIMING, 'TIMING')
        log.addLevelName(self.RUN, 'RUN')
        log.basicConfig(format='[%(levelname)s] %(message)s')

    def debug(self, msg, *args, **kwargs):
        if self.isEnabledFor(self.DEBUG):
            self._log(self.DEBUG, msg, args, **kwargs)

    def info(self, msg, *args, **kwargs):
        if self.isEnabledFor(self.INFO) or self.isEnabledFor(self.RUN):
            self._log(self.INFO, msg, args, **kwargs)

    def warning(self, msg, *args, **kwargs):
        print(self.isEnabledFor(self.RUN))
        if self.isEnabledFor(self.WARNING) and self.isEnabledFor(self.INFO) or not self.isEnabledFor(self.RUN):
            self._log(self.WARNING, msg, args, **kwargs)

    def timing(self, msg, *args, **kwargs):
        if self.isEnabledFor(self.TIMING):
            self._log(self.TIMING, msg, args, **kwargs)

log.setLoggerClass(myLogger)
logging = log.getLogger("main")








