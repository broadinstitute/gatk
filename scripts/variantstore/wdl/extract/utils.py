import time
import google.api_core.exceptions
from google.cloud import bigquery


def execute_with_retry(client, label, sql):
    """Run a BigQuery SQL string with a label.

    Three retries with incremental backoff if BiqQuery returns 'retry-able errors'
    (see https://googleapis.dev/python/bigquery/latest/_modules/google/api_core/retry.html),
    any other error is re-raised.

    Parameters
    ----------
    client : bigquery.Client
        with credentials, project, and default_query_job_config already defined
    label : str
        additional label to add to job
    sql : str
        SQL to run
    """
    retry_delay = [30, 60, 90]
    start = time.time()
    while len(retry_delay) >= 0:
        try:
            query_label = label.replace(" ", "-").strip().lower()
            existing_labels = client._default_query_job_config.labels
            job_labels = existing_labels
            job_labels["gvs_query_name"] = query_label
            job_config = bigquery.QueryJobConfig(labels=job_labels)
            query = client.query(sql, job_config=job_config)
            print(f"STARTING - {label} (jobid: {query.job_id})")
            results = query.result()
            job = client.get_job(query.job_id)
            mb_billed = int(0 if job.total_bytes_billed is None else job.total_bytes_billed) / (
                        1024 * 1024)
            print(
                f"COMPLETED ({time.time() - start} seconds, {3 - len(retry_delay)} retries, {mb_billed} MBs) - {label}")

            return {'results': results, 'job': job, 'label': label}
        except (google.api_core.exceptions.InternalServerError,
                google.api_core.exceptions.TooManyRequests,
                google.api_core.exceptions.ServiceUnavailable) as err:
            if len(retry_delay) > 0:
                t = retry_delay.pop(0)
                print(f"Error {err} running query {label}, sleeping for {t}")
                time.sleep(t)
            else:
                raise err
        except Exception:
            raise


def utf8len(s):
    return len(s.encode('utf-8'))


def write_job_stats(jobs, client, fq_dataset, call_set_identifier, step, call, shard_identifier):
    total = 0

    for query_return in jobs:
        job = client.get_job(query_return['job'])
        bytes_billed = int(0 if job.total_bytes_billed is None else job.total_bytes_billed)
        total = total + bytes_billed
        print(query_return['label'], "-- jobid:  (", job.job_id, ") <====> Cache Hit:", job.cache_hit, bytes_billed/(1024 * 1024), " MiBs")
    print("\nTotal GiBs billed ", total/(1024 * 1024 * 1024), " GiBs\n")

    # populate cost_observability data
    sql = f"""INSERT INTO `{fq_dataset}.cost_observability`
            (call_set_identifier, step, call, shard_identifier, event_key, call_start_timestamp, event_timestamp, event_bytes)
            VALUES('{call_set_identifier}', '{step}', '{call}', '{shard_identifier}', 'BigQuery Query Scanned',
            CURRENT_TIMESTAMP(), CURRENT_TIMESTAMP(), {total})"""
    query = client.query(sql)
    query.result()
