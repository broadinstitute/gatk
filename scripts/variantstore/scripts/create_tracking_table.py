#!/usr/bin/env python3
"""
DEPRECATED: The parquet_load_status tracking table has been eliminated.

Idempotency is now determined by inspecting BigQuery directly
(INFORMATION_SCHEMA.PARTITIONS for vet/ref_ranges tables, and the
sample_chromosome_ploidy table itself for ploidy data) rather than
maintaining a separate tracking table.

This script is retained as a no-op stub for backward compatibility and
will be removed in a future cleanup.
"""

import sys


def main():
    print(
        "NOTE: create_tracking_table.py is deprecated and has no effect. "
        "The parquet_load_status table is no longer used by the GVS Parquet "
        "loading pipeline. Idempotency is now determined by querying BigQuery "
        "directly via INFORMATION_SCHEMA.PARTITIONS."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
