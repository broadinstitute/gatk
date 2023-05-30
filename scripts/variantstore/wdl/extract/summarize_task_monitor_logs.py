import argparse
import sys
import re
import os

MaxCpu = -100.0
MaxMem = -100.0
MaxMemPct = -100.0
MaxDisk = -100.0
MaxDiskPct = -100.0


def parse_monitoring_log_files(mlog_files, output_file):
    with open(output_file, 'w') as output:
        header = f"Total Mem\tMax Mem Used\tMax Mem Used (%)\tTotal Disk\tMax Disk Used\tMax Disk Used (" \
                 f"%)\tTask\tShard\tFile\n"
        output.write(header)

        for mlog_file in mlog_files:
            if not os.path.exists(mlog_file):
                eprint(f"ERROR: File {mlog_file} does not exist")
            parse_monitoring_log_file(mlog_file, output)


def parse_monitoring_log_file(mlog_file, output):
    eprint(f"Parsing: {mlog_file}")

    if os.stat(mlog_file).st_size == 0:
        eprint(f"Skipping zero length file")
        return

    with open(mlog_file) as ml:
        advance_to_section(ml, "--- General Information")

        line = ml.readline().rstrip()  # #CPU: 16
        tokens = line.split(": ")
        if tokens[0] != '#CPU':
            eprint(f"ERROR: Line '{line}' does not look as expected. Is this a monitoring_log file?")
            sys.exit(1)
        num_cpus = tokens[1]
        eprint(f"Num CPUs: {num_cpus}")

        line = ml.readline().rstrip()  # Total Memory: 98.25 GiB
        tokens = line.split(": ")
        if tokens[0] != 'Total Memory':
            eprint(f"ERROR: Line '{line}' does not look as expected. Is this a monitoring_log file?")
            sys.exit(1)
        subtokens = tokens[1].split()
        if len(subtokens) != 2:
            eprint(f"ERROR: Line '{line}' does not look as expected. Is this a monitoring_log file?")
            sys.exit(1)
        total_memory = float(subtokens[0])
        total_memory_units = subtokens[1]
        eprint(f"Total Memory: {total_memory} {total_memory_units}")

        line = ml.readline().rstrip()  # Total Disk space: 985.000 GiB
        tokens = line.split(": ")
        if tokens[0] != 'Total Disk space':
            eprint(f"ERROR: Line '{line}' does not look as expected. Is this a monitoring_log file?")
            sys.exit(1)
        subtokens = tokens[1].split()
        if len(subtokens) != 2:
            eprint(f"ERROR: Line '{line}' does not look as expected. Is this a monitoring_log file?")
            sys.exit(1)
        total_disk = float(subtokens[0])
        total_disk_units = subtokens[1]
        eprint(f"Total Disk space: {total_disk} {total_disk_units}")

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

        field_type = 0
        while line := ml.readline().rstrip():
            if field_type == 1:
                parse_cpu_usage_line(line)
            elif field_type == 2:
                parse_memory_usage_line(line)
            elif field_type == 3:
                parse_disk_usage_line(line)
            field_type += 1
            if field_type > 4:
                field_type = 0

        file_path = os.path.abspath(mlog_file)
        tokens = file_path.split("/")
        # Looks like '*/call-IndelsVariantRecalibrator/monitoring.log' for non-sharded
        # Looks like '*/call-ScoreVariantAnnotationsINDELs/shard-35/monitoring.log' for sharded
        # Looks like '*/call-MergeVCFs/cacheCopy/monitoring.log' for cached, non-sharded
        # Looks like '*/call-ExtractFilterTask/shard-0/cacheCopy/monitoring.log' for cached sharded
        shard = ""
        index = -2
        if tokens[index] == "cacheCopy":
            index = -3
        if tokens[index].startswith("shard-"):
            shard = tokens[index][6:]
            task = tokens[index - 1][5:]  # Strip off the 'call-' prefix
        else:
            task = tokens[index][5:]

        summary = f"{total_memory}\t{MaxMem}\t{MaxMemPct}\t{total_disk}\t{MaxDisk}\t{MaxDiskPct}\t{task}\t{shard}\t" \
                  f"{os.path.abspath(mlog_file)}\n"
        output.write(summary)


def parse_cpu_usage_line(line):
    p = "^\\* CPU usage\\: (\\d+\\.\\d+)%$"  # Looks Like: * CPU usage: 17.6%
    m = re.match(p, line)
    if m is not None:
        cpu = float(m.group(1))
        global MaxCpu
        if cpu > MaxCpu:
            MaxCpu = cpu
    else:
        # Check if it's a nan (we see this sometimes at startup)
        p2 = "^\\* CPU usage\\: -nan\\%$"  # * CPU usage: -nan%
        m2 = re.match(p2, line)
        if m2 is None:
            # Check if it's just empty
            p3 = "^\\* CPU usage\\: *$"  # * CPU usage:
            m3 = re.match(p3, line)
            if m3 is None:
                eprint(f"ERROR: Line '{line}' does not look like a CPU usage line. Is this a monitoring_log file?")
                sys.exit(1)


def parse_memory_usage_line(line):
    p = "^\\* Memory usage\\: (\\d+\\.\\d+) \\S+ (\\d+(\\.\\d*)?)\\%$"  # Looks Like: * Memory usage: 1.79 GiB 1.8%
    m = re.match(p, line)
    if m is None:
        eprint(f"ERROR: Line '{line}' does not look like a Memory usage line. Is this a monitoring_log file?")
        sys.exit(1)
    mem = float(m.group(1))
    mem_pct = float(m.group(2))
    global MaxMem
    global MaxMemPct
    if mem > MaxMem:
        MaxMem = mem
        MaxMemPct = mem_pct


def parse_disk_usage_line(line):
    p = "^\\* Disk usage\\: (\\d+\\.\\d+) \\S+ (\\d+(\\.\\d*)?)\\%$"  # Looks Like: * Disk usage: 22.000 GiB 3%
    m = re.match(p, line)
    if m is None:
        eprint(f"ERROR: Line '{line}' does not look like a Disk usage line. Is this a monitoring_log file?")
        sys.exit(1)
    disk = float(m.group(1))
    disk_pct = float(m.group(2))
    global MaxDisk
    global MaxDiskPct
    if disk > MaxDisk:
        MaxDisk = disk
        MaxDiskPct = disk_pct


def advance_to_section(fp, section_header_start):
    """
    A method to advance to the runtime block in the file
    :param fp: The file's file pointer
    :param section_header_start: The (starting part) of the section header in the file.
    :return: The nextmost line that is not a header line
    """
    while line := fp.readline():
        if line.startswith(section_header_start):
            return line


def eprint(*the_args, **kwargs):
    print(*the_args, file=sys.stderr, **kwargs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='A tool to summarize the output of multiple monitoring logs')

    parser.add_argument('--input', nargs='+', help='Monitoring log file(s)', required=True)
    parser.add_argument('--output', type=str, help='Output Monitoring log summary file', required=True)
    args = parser.parse_args()

    parse_monitoring_log_files(args.input, args.output)
