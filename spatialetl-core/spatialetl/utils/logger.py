# MIT License
# Copyright (c) 2024 [SNALE - French SAS Company - RCS 951 724 616]
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
import logging as log
class RunFilter(log.Filter):
    def filter(self, record):
        if record.levelno==log.WARN or record.levelno==myLogger.TIMING:
            return False
        return True

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
        if self.isEnabledFor(self.WARNING):
            self._log(self.WARNING, msg, args, **kwargs)

    def timing(self, msg, *args, **kwargs):
        if self.isEnabledFor(self.TIMING):
            self._log(self.TIMING, msg, args, **kwargs)

    def setLevel(self,level):
        log.getLogger().setLevel(level)
        if level == self.RUN:
            self.addFilter(RunFilter())

log.setLoggerClass(myLogger)
logging = log.getLogger("main")









