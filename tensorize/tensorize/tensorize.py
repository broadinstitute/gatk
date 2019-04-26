import apache_beam as beam
from apache_beam.options.pipeline_options import PipelineOptions

from .utils import count_ones, process_entries, write_tensor_from_sql


def run(pipeline_options: PipelineOptions, output_file: str):
    p = beam.Pipeline(options=pipeline_options)

    bigquery_source = beam.io.BigQuerySource(
        query='select * from `lubitz.coding`',
        use_standard_sql=True
    )

    # Query table in BQ
    table_data = (
        p
        | 'QueryTable' >> beam.io.Read(bigquery_source)
        # Each row is a dictionary where the keys are the BigQuery columns

        | 'CreateKey' >> beam.Map(lambda elem: (elem['coding_file_id'], 1))
        # group by key

        | 'GroupByKey' >> beam.GroupByKey()
        # count entries per key

        | 'Count' >> beam.Map(count_ones)
    )

    # Should be replaced by the schema.json
    # table_schema = bigquery.TableSchema()
    # coding_file_id_field = bigquery.TableFieldSchema()
    # coding_file_id_field.name = 'coding_file_id'
    # coding_file_id_field.type = 'INTEGER'
    # table_schema.fields.append(coding_file_id_field)
    #
    # coding_field = bigquery.TableFieldSchema()
    # coding_field.name = 'coding'
    # coding_field.type = 'STRING'
    # table_schema.fields.append(coding_field)
    #
    # meaning_field = bigquery.TableFieldSchema()
    # meaning_field.name = 'meaning'
    # meaning_field.type = 'STRING'
    # table_schema.fields.append(meaning_field)

    # write_table_spec = bigquery.TableReference(
    #     projectId='broad-ml4cvd',
    #     datasetId='pb_working',
    #     tableId=RUN_NAME
    # )
    #
    # write transformed to new table
    # table_data | 'Writing to BigQuery' >> beam.io.WriteToBigQuery(
    #    write_table_spec,
    #    schema=table_schema,
    #    write_disposition=beam.io.BigQueryDisposition.WRITE_TRUNCATE,
    #    create_disposition=beam.io.BigQueryDisposition.CREATE_IF_NEEDED
    # )

    # Write to file
    table_data | 'WriteToFile' >> beam.io.WriteToText(output_file)

    result = p.run()
    result.wait_until_finish()


def run2(pipeline_options: PipelineOptions, output_file: str):
    limit = 500

    p = beam.Pipeline(options=pipeline_options)

    bigquery_source = beam.io.BigQuerySource(
        # query='select * from `ukbb7089_r10data.phenotype` limit %s' % limit,
        query='select * from `ukbb7089_r10data.phenotype`',
        use_standard_sql=True
    )

    # Query table in BQ
    table_data = (
            p
            | 'QueryTable' >> beam.io.Read(bigquery_source)

            # Each row is a dictionary where the keys are the BigQuery columns
            | 'CreateKey' >> beam.Map(lambda row: (row['sample_id'], row))

            # Group by key
            | 'GroupByKey' >> beam.GroupByKey()

            # Create dict of sample_id -> list(row dicts)
            | 'ProcessSamples' >> beam.Map(process_entries)
    )

    # Write to file
    table_data | 'WriteToFile' >> beam.io.WriteToText(output_file)

    result = p.run()
    result.wait_until_finish()


def tensorize_categorical_fields(pipeline_options: PipelineOptions):
    p = beam.Pipeline(options=pipeline_options)

    # TODO: Don't hardcode LIMIT in the query
    limit = 300
    query = """
        SELECT d.field, p_sub_f_c.meaning, p_sub_f_c.sample_id, p_sub_f_c.fieldid, p_sub_f_c.instance, p_sub_f_c.array_idx, p_sub_f_c.value FROM
        (
            SELECT c.meaning, p_sub.sample_id, p_sub.fieldid, p_sub.instance, p_sub.array_idx, p_sub.value FROM
                (
                    SELECT p.sample_id, p.fieldid, p.instance, p.array_idx, p.value, p.coding_file_id FROM `ukbb_dev.phenotype` p
                    INNER JOIN (SELECT * FROM `shared_data.tensorization_fieldids` LIMIT 300) f
                    ON p.fieldid = f.fieldid
                ) AS p_sub
            JOIN `ukbb_dev.coding` c
            ON c.coding=p_sub.value AND c.coding_file_id=p_sub.coding_file_id
        ) AS p_sub_f_c
        JOIN `ukbb_dev.dictionary` d
        ON p_sub_f_c.fieldid = d.fieldid"""
    
    bigquery_source = beam.io.BigQuerySource(query=query, use_standard_sql=True)

    # Query table in BQ
    steps = (
            p
            | 'QueryTables' >> beam.io.Read(bigquery_source)

            # Each row is a dictionary where the keys are the BigQuery columns
            | 'CreateKey' >> beam.Map(lambda row: (row['sample_id'], row))

            # Group by key
            | 'GroupByKey' >> beam.GroupByKey()

            # Create dict of sample_id -> list(row dicts)
            | 'ProcessSamples' >> beam.Map(process_entries)

            # Format into hd5 files and upload to GCS
            | 'CreateHd5sAndUploadToGCS' >> beam.Map(write_tensor_from_sql)
    )

    result = p.run()
    result.wait_until_finish()
