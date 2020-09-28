"""Methods for storing and retrieving annotations within notebooks."""

import abc
import datetime
from typing import Optional, Union

from google.cloud import bigquery
from google.cloud.bigquery import magics as bqmagics
import pandas as pd


class AnnotationStorage(abc.ABC):
  """Base class for annotation storage.

  This strategy pattern to allow different storage mechanisms to be used in different notebooks environments.
  """

  @abc.abstractmethod
  def describe(self) -> str:
    """Return a string describing how annotations are stored."""

  @abc.abstractmethod
  def submit_annotation(
      self, sample_id: Union[int, str], annotator: str, key: str,
      value_numeric: Optional[Union[int, float]], value_string: Optional[str], comment: str,
  ) -> bool:
    """Add an annotation to the collection of annotations.

    Args:
      sample_id: Id of the sample this annotation concerns.
      annotator: Email address or hostname associated with the person submitting the annotation.
      key: The key of the annotation, such as a phenotype field name.
      value_numeric: The numeric value of the field being annotated, if applicable.
      value_string: The string value of the field being annotated, if applicable.
      comment: The actual annotation.
    Returns:
      Whether the submission was successful. Throws an Exception on failure.
    """

  @abc.abstractmethod
  def view_recent_submissions(self, count: int = 10) -> pd.DataFrame:
    """View a dataframe of up to [count] most recent submissions.

    Args:
      count: The number of the most recent submissions to return.

    Returns:
      A dataframe of the most recent annotations.
    """


class TransientAnnotationStorage(AnnotationStorage):
  """Store annotations temporarily in memory.

  This is useful for demonstration purposes.
  """

  def __init__(self):
    self.annotations = []

  def describe(self) -> str:
    return '''Annotations will be stored in memory only during the duration of this demo.\n
    For durable storage of annotations, use BigQueryAnnotationStorage instead.'''

  def submit_annotation(
      self, sample_id: Union[int, str], annotator: str, key: str,
      value_numeric: Optional[Union[int, float]], value_string: Optional[str], comment: str,
  ) -> bool:
    """Add this annotation to our in-memory collection of annotations.

    Args:
      sample_id: Id of the sample this annotation concerns.
      annotator: Email address or hostname associated with the person submitting the annotation.
      key: The key of the annotation, such as a phenotype field name.
      value_numeric: The numeric value of the field being annotated, if applicable.
      value_string: The string value of the field being annotated, if applicable.
      comment: The actual annotation.
    Returns:
      True
    """
    annotation = {
        'sample_id': sample_id,
        'annotator': annotator,
        'annotation_timestamp': datetime.datetime.now(),
        'key': key,
        'value_numeric': value_numeric,
        'value_string': value_string,
        'comment': comment,
    }
    self.annotations.append(annotation)
    return True

  def view_recent_submissions(self, count: int = 10) -> pd.DataFrame:
    """View a dataframe of up to [count] most recent submissions.

    Args:
      count: The number of the most recent submissions to return.

    Returns:
      A dataframe of the most recent annotations.
    """
    return pd.DataFrame.from_dict(self.annotations[-1 * count :])


class BigQueryAnnotationStorage(AnnotationStorage):
  """Store annotations in a BigQuery table.

  The table must have the schema in file annotations_schema.json.
  For example, the table can be created with the bq command line tool:

    bq --project your-project-id mk \\
      --table \\
      --description 'Annotations table for Terra ml4h featured workspace' \\
      your-dataset.annotations \\
      annotations_schema.json
  """

  def __init__(self, table: str):
    """This table should already exist."""
    self.table = table

  def describe(self) -> str:
    return f'''Annotations are stored in BigQuery table {self.table}'''

  def submit_annotation(
      self, sample_id: Union[int, str], annotator: str, key: str,
      value_numeric: Optional[Union[int, float]], value_string: Optional[str], comment: str,
  ) -> bool:
    """Call a BigQuery INSERT statement to add a row containing annotation information.

    Args:
      sample_id: Id of the sample this annotation concerns.
      annotator: Email address or hostname associated with the person submitting the annotation.
      key: The key of the annotation, such as a phenotype field name.
      value_numeric: The numeric value of the field being annotated, if applicable.
      value_string: The string value of the field being annotated, if applicable.
      comment: The actual annotation.
    Returns:
      Whether the submission is complete. Throws an Exception on failure.
    """
    # Set up biquery client.
    bqclient = bigquery.Client(credentials=bqmagics.context.credentials)

    # Format the insert string.
    query_string = f'''
        INSERT INTO `{self.table}`
          (sample_id, annotator, annotation_timestamp, key, value_numeric, value_string, comment)
        VALUES
          ('{sample_id}', '{annotator}', CURRENT_TIMESTAMP(), '{key}', SAFE_CAST('{value_numeric}' as NUMERIC), '{value_string}', '{comment}')
        '''

    # Submit the insert request.
    submission = bqclient.query(query_string)

    # Check for any errors. Upon error, this will throw an exception.
    _ = submission.result()

    # Return whether the submission completed.
    return submission.done()

  def view_recent_submissions(self, count: int = 10) -> pd.DataFrame:
    """View a dataframe of up to [count] most recent submissions.

    This is a convenience method for use within the annotation flow. For full access to the underlying annotations,
    connect to the table directly using the BigQuery client of your choice.

    Args:
      count: The number of the most recent submissions to return.

    Returns:
      A dataframe of the most recent annotations.
    """

    # set up biquery client.
    bqclient = bigquery.Client(credentials=bqmagics.context.credentials)

    # Format the query string.
    query_string = f'''
        SELECT * FROM `{self.table}`
        ORDER BY annotation_timestamp DESC
        LIMIT {count}
        '''

    # submit the query and store the result as a dataframe
    df = bqclient.query(query_string).result().to_dataframe()

    return df
