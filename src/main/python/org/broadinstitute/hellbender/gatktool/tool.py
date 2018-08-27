"""
Functions and classes used to extend a GATK tool with Python.

GATK uses two FIFOs to communicate wth Python. The "ack" FIFO is read by GATK
and written by Python code, and is used to signal that a Python command has
completed execution. The "data" FIFO is written by GATK and read by Python,
and is used to pass data to Python from Java.

Most of the functions in this module are intended to be called by GATK via
the StreamingPythonScriptExecutor Java class, and are not called by Python
code directly. The one exception is the readDataFIFO function, which can be
used to read data that had been passed to Python by GATK Java code.
"""

import sys
import os
import cProfile, pstats, io

_ackFIFO = None
_dataFIFO = None
_GATKProfiler = None

def initializeGATK(ackFIFOName: str):
    """
    Open the GATK ack FIFO and install the exception handler hook.

    Called by GATK when the StreamingPythonScriptExecutor is initialized,
    which is normally in onTraversalStart. Initializes the ack FIFO and
    installs the exception hook. Since the exception hook uses the ack FIFO,
    it can't be installed until after the FIFO is initialized.
    """
    global _ackFIFO
    _ackFIFO = AckFIFO(ackFIFOName)
    sys.excepthook = gatkExceptionHook


def gatkExceptionHook(exceptionType, value, traceback):
    """
    GATK Handler for uncaught Python exceptions.

    The is installed by initializeGATK after the ack FIFO has been
    initialized. When an unhandled exception is caught, the handler
    sends a nack to GATK through the FIFO, which results in a
    PythonScriptExecutorException being thrown in the tool.
    """
    sendNack()
    sys.__excepthook__(exceptionType, value, traceback)


def sendAck():
    """
    Send a positive acknowledgment to GATK. This should generally only
    be called by python code that is embedded in Java, since the executor
    keeps track of whether an ack request is outstanding.
    """
    global _ackFIFO
    _ackFIFO.writeAck()


def sendNack():
    """
    Send a negative acknowledgment to GATK. Generally only called by the
    installed exception hook. This will result in a Java exception being
    thrown that unless caught by Java code, will terminate the tool.
    """
    global _ackFIFO
    _ackFIFO.writeNack()


def terminateGATK():
    """
    Called by GATK when no more Python commands will be executed
    """
    global _ackFIFO
    if _ackFIFO is None:
        raise RuntimeError("ack FIFO has not been initialized")
    _ackFIFO.close()
    _ackFIFO = None


def initializeDataFIFO(dataFIFOName: str):
    """
    Initialize the data FIFO for reading.

    Once this method has been called, the FIFO may be read using the
    readDataFIFO function.
    """
    global _dataFIFO
    _dataFIFO = DataFIFO(dataFIFOName)


def closeDataFIFO():
    """"
    Close the data FIFO. Afteer this method hass been called, the
    data FIFO can no longer be read.
    """
    global _dataFIFO
    if _dataFIFO is None:
        raise RuntimeError("data FIFO has not been initialized")
    _dataFIFO.close()
    _dataFIFO = None


def readDataFIFO() -> str:
    """
    Read a line from the Data FIFO.
    :return: string
    """
    global _dataFIFO
    return _dataFIFO.readLine()

def startProfiling():
    """
    Start Python CProfile profiling.
    """
    global _GATKProfiler
    _GATKProfiler = cProfile.Profile()
    _GATKProfiler.enable()

def endProfiling(profileName: str):
    """
    End Python CProfile profiling and write results to a file. The
    startProfile function must have been previously called. The results
    are ordered by cummulative time.
    :param profileName: name of the file to which the profiling results should be written.
    """
    global _GATKProfiler
    _GATKProfiler.disable()
    gatkProfilerDescriptor = os.open(profileName, os.O_WRONLY | os.O_CREAT)
    gatkProfileStream = os.fdopen(gatkProfilerDescriptor, 'w')
    gatkStats = pstats.Stats(_GATKProfiler, stream=gatkProfileStream).sort_stats('cumulative')
    gatkStats.print_stats()
    gatkProfileStream.close()
    del gatkProfileStream
    del gatkProfilerDescriptor
    del gatkStats


class AckFIFO:
    """
    Manage the FIFO used to notify GATK (via an ack) that a command has
    completed, or failed due to an unhandled exception (via a nck).
    """
    _ackString = "ack"
    _nackString = "nck"

    def __init__(self, ackFIFOName: str) -> None:
        """Open the ack fifo stream for writing only"""
        self.ackFIFOName = ackFIFOName
        writeDescriptor = os.open(self.ackFIFOName, os.O_WRONLY)
        self.fileWriter = os.fdopen(writeDescriptor, 'w')

    def writeAck(self):
        """
        Write a positive acknowledgement to the ACK FIFO.
        """
        if self.fileWriter is None:
            raise RuntimeError("ack FIFO has not been initialized")
        self.fileWriter.write(AckFIFO._ackString)
        self.fileWriter.flush()

    def writeNack(self):
        """
        Write a negative acknowledgement to the ACK FIFO.

        Calling this method will result in an exception being thrown
        in the GATK tool on whos behalf this module is running.
        """
        if self.fileWriter is None:
            raise RuntimeError("ack FIFO has not been initialized")
        self.fileWriter.write(AckFIFO._nackString)
        self.fileWriter.flush()

    def close(self):
        assert self.fileWriter != None
        self.fileWriter.close()
        self.fileWriter = None


class DataFIFO:
    """
    Manage the FIFO stream used for transferring data from the GATK tool to
    Python code.

    The FIFO is written by GATK and read by Python.
    """
    def __init__(self, dataFIFOName: str) -> None:
        """Open the data stream fifo for reading"""
        self.dataFIFOName = dataFIFOName

        # the data fifo is always opened for read only on the python side
        readDescriptor = os.open(self.dataFIFOName, os.O_RDONLY)
        self.fileReader = os.fdopen(readDescriptor, 'r')

    def readLine(self) -> str:
        """
        Read a single line from the Data FIFO.
        :return: string
        """
        if self.fileReader is None:
            raise RuntimeError("data FIFO reader has not been initialized")
        return self.fileReader.readline()

    def close(self):
        if self.fileReader is None:
            raise RuntimeError("data FIFO reader has not been initialized")
        self.fileReader.close()
        self.fileReader = None
