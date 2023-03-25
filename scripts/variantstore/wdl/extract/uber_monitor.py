import argparse
import sys
import re
import os

global MaxCpu
global MaxMem
global MaxMemPct
global MaxDisk
global MaxDiskPct

def parse_monitoring_log_files(mlog_files):
    header = f"Total Mem\tMax Mem Used\tMax Mem Used (%)\tTotal Disk\tMax Disk Used\tMax Disk Used (%)\tFile"
    print(header)

    for mlog_file in mlog_files:
        parse_monitoring_log_file(mlog_file)

def parse_monitoring_log_file(mlog_file):
    if (args.verbose):
        eprint(f"Parsing: {mlog_file}")

    if (os.stat(mlog_file).st_size == 0):
        if (args.verbose):
            eprint(f"Skipping zero length file")
        return

    with open(mlog_file) as ml:
        advance_to_section(ml, "--- General Information")

        line = ml.readline().rstrip()       # #CPU: 16
        tokens = line.split(": ")
        if (tokens[0] != '#CPU'):
            eprint(f"ERROR: Line '{line}' does not look as expected. Is this a monitoring_log file?")
            sys.exit(1)
        NumCPUS = tokens[1]
        if (args.verbose):
            eprint(f"Num CPUs: {NumCPUS}")

        line = ml.readline().rstrip()       # Total Memory: 98.25 GiB
        tokens = line.split(": ")
        if (tokens[0] != 'Total Memory'):
            eprint(f"ERROR: Line '{line}' does not look as expected. Is this a monitoring_log file?")
            sys.exit(1)
        subtokens = tokens[1].split()
        if (len(subtokens) != 2):
            eprint(f"ERROR: Line '{line}' does not look as expected. Is this a monitoring_log file?")
            sys.exit(1)
        TotalMemory = float(subtokens[0])
        TotalMemoryUnits = subtokens[1]
        if (args.verbose):
            eprint(f"Total Memory: {TotalMemory} {TotalMemoryUnits}")

        line = ml.readline().rstrip()       # Total Disk space: 985.000 GiB
        tokens = line.split(": ")
        if (tokens[0] != 'Total Disk space'):
            eprint(f"ERROR: Line '{line}' does not look as expected. Is this a monitoring_log file?")
            sys.exit(1)
        subtokens = tokens[1].split()
        if (len(subtokens) != 2):
            eprint(f"ERROR: Line '{line}' does not look as expected. Is this a monitoring_log file?")
            sys.exit(1)
        TotalDisk = float(subtokens[0])
        TotalDiskUnits = subtokens[1]
        if (args.verbose):
            eprint(f"Total Disk space: {TotalDisk} {TotalDiskUnits}")

        advance_to_section(ml, "--- Runtime Information")

        global MaxCpu
        MaxCpu = -100.0
        global MaxMem
        MaxMem = -100.0
        global MaxMemPct
        MaxMemPct = -100.0
        global MaxDisk
        MaxDisk = -100.0
        global MaxDiskPct
        MaxDiskPct = -100.0

        type = 0
        while (line := ml.readline().rstrip()):
            if type == 0:
                parse_timestamp_line(line)
            elif type == 1:
                parse_cpu_usage_line(line)
            elif type == 2:
                parse_memory_usage_line(line)
            elif type == 3:
                parse_disk_usage_line(line)
            elif type == 4:
                parse_rwio_line(line)
            type += 1
            if (type > 4):
                type = 0

        summary = f"{TotalMemory}\t{MaxMem}\t{MaxMemPct}\t{TotalDisk}\t{MaxDisk}\t{MaxDiskPct}\t{os.path.abspath(mlog_file)}"
        print(summary)

def parse_timestamp_line(line):
    line = ""
    # Do nothing for now.   ## [Tue Mar 14 14:52:39 UTC 2023]

def parse_cpu_usage_line(line):
    p = "^\* CPU usage\: (\d+\.\d+)\%$"     #* CPU usage: 17.6%
    m = re.match(p, line)
    if (m != None):
        cpu = float(m.group(1))
        if (args.verbose):
            eprint(f"CPU: {cpu}")
        global MaxCpu
        if (cpu > MaxCpu):
            MaxCpu = cpu
    else:
        # Check if it's a nan (see this sometimes at startup)
        p2 = "^\* CPU usage\: -nan\%$"     #* CPU usage: -nan%
        m2 = re.match(p2, line)
        if (m2 == None):
            eprint(f"ERROR: Line '{line}' does not look like a CPU usage line. Is this a monitoring_log file?")
            sys.exit(1)

def parse_memory_usage_line(line):
    p = "^\* Memory usage\: (\d+\.\d+) \S+ (\d+(\.\d*)?)\%$"        #* Memory usage: 1.79 GiB 1.8%
    m = re.match(p, line)
    if (m == None):
        eprint(f"ERROR: Line '{line}' does not look like a Memory usage line. Is this a monitoring_log file?")
        sys.exit(1)
    if (args.verbose):
        eprint(f"Memory: {m.group(1)} {m.group(2)}")
    mem = float(m.group(1))
    memPct = float(m.group(2))
    global MaxMem
    global MaxMemPct
    if (mem > MaxMem):
        MaxMem = mem
        MaxMemPct = memPct

def parse_disk_usage_line(line):
    p = "^\* Disk usage\: (\d+\.\d+) \S+ (\d+(\.\d*)?)\%$"      #* Disk usage: 22.000 GiB 3%
    m = re.match(p, line)
    if (m == None):
        eprint(f"ERROR: Line '{line}' does not look like a Disk usage line. Is this a monitoring_log file?")
        sys.exit(1)
    if (args.verbose):
        eprint(f"Disk: {m.group(1)} {m.group(2)}")
    disk = float(m.group(1))
    diskPct = float(m.group(2))
    global MaxDisk
    global MaxDiskPct
    if (disk > MaxDisk):
        MaxDisk = disk
        MaxDiskPct = diskPct

def parse_rwio_line(line):
    line = ""
    # Do nothing for now.       # * Read/Write IO: 0.000 MiB/s 127.114 MiB/s

def advance_to_section(fp, section_header_start):
    """
    A method to advance to the runtime block in the file
    :param fp: The file's file pointer
    :return: The nextmost line that is not a header line
    """
    while (line := fp.readline()):
        if line.startswith(section_header_start):
            return line

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='A tool to parse the output of multiple monitoring logs')

    parser.add_argument('--input', nargs='+', help='Monitoring log file(s)', required=True)
    parser.add_argument("--verbose", action="store_true", help="increase output verbosity")
    args = parser.parse_args()

    parse_monitoring_log_files(args.input)