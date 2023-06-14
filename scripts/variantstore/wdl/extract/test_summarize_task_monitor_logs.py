import unittest

from summarize_task_monitor_logs import parse_monitoring_log_files


class TestSummarizeTaskMonitorLogs(unittest.TestCase):
    def test_parse_monitoring_log_files(self):
        import tempfile

        input_files = ['summarize_task_monitor_logs_test_files/call-ScoreVariantAnnotationsINDELs/shard-34/monitoring'
                       '.log',
                       'summarize_task_monitor_logs_test_files/call-ScoreVariantAnnotationsINDELs/shard-35/monitoring'
                       '.log',
                       'summarize_task_monitor_logs_test_files/call-IndelsVariantRecalibrator/monitoring.log',
                       'summarize_task_monitor_logs_test_files/call-ExtractFilterTask/shard-0/cacheCopy/monitoring.log',
                       'summarize_task_monitor_logs_test_files/call-MergeVCFs/cacheCopy/monitoring.log',
                       'summarize_task_monitor_logs_test_files/call-SamplesTableDatetimeCheck/monitoring.log']
        with tempfile.NamedTemporaryFile() as actual_output_file:
            parse_monitoring_log_files(input_files, actual_output_file.name)

            expected_output_file = 'summarize_task_monitor_logs_test_files/expected_monitoring_summary_file.txt'
            with open(actual_output_file.name, 'r') as actual, open(expected_output_file, 'r') as expected:
                while True:
                    actual_line = actual.readline().rstrip()
                    expected_line = expected.readline().rstrip()
                    if (actual_line == "") or (expected_line == ""):
                        break
                    # We Need to ignore any diffs in the file path (last element)
                    # as it differs on the test running under the docker.
                    actual_line = ",".join(actual_line.split("\t")[0:-1])
                    expected_line = ",".join(expected_line.split("\t")[0:-1])
                    self.assertEqual(actual_line, expected_line)
