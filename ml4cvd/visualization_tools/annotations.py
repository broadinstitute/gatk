"""Methods for capturing and displaying annotations within notebooks."""

import os
import socket
from IPython.display import display
from IPython.display import HTML
import ipywidgets as widgets
from ml4cvd.visualization_tools.annotation_storage import AnnotationStorage
from ml4cvd.visualization_tools.annotation_storage import TransientAnnotationStorage

DEFAULT_ANNOTATION_STORAGE = TransientAnnotationStorage()


def _get_df_sample(sample_info, sample_id):
  """Return a dataframe containing only the row for the indicated sample_id."""
  df_sample = sample_info[sample_info['sample_id'] == str(sample_id)]
  if 0 == df_sample.shape[0]: df_sample = sample_info.query('sample_id == ' + str(sample_id))
  return df_sample


def display_annotation_collector(sample_info, sample_id, annotation_storage: AnnotationStorage = DEFAULT_ANNOTATION_STORAGE, custom_annotation_key=None):
  """Method to create a gui (set of widgets) through which the user can create an annotation and submit it to storage.

  Args:
    sample_info: Dataframe containing tabular data for all the samples.
    sample_id: The selected sample for which the values will be displayed.
    annotation_storage: An instance of AnnotationStorage.
    custom_annotation_key: The key for an annotation of data other than the tabular fields.

  Returns:
    A notebook-friendly messages indicating the status of the submission.
  """

  df_sample = _get_df_sample(sample_info, sample_id)
  if df_sample.shape[0] == 0:
    return HTML(f'''<div class="alert alert-block alert-danger">
      <b>Warning:</b> Sample {sample_id} not present in sample_info DataFrame.</div>''')

  # Show the sample ID for this annotation.
  sample = widgets.HTML(value=f'For sample <b>{sample_id}</b>')

  # Allow the user to pick a key about which to comment.
  annotation_keys = []
  if custom_annotation_key:
    annotation_keys.append(custom_annotation_key)
  annotation_keys.extend(sorted(sample_info.keys()))
  key = widgets.Dropdown(
      options=annotation_keys,
      description='Key:',
      disabled=False,
  )

  # Return the sample's value for that key, when applicable.
  valuelabel = widgets.Label(value='Value: ')
  keyvalue = widgets.Label(value=None if custom_annotation_key else str(df_sample[key.value].iloc[0]))

  box1 = widgets.HBox(
      [key, valuelabel, keyvalue],
      layout=widgets.Layout(width='50%'),
  )

  # Have keyvalue auto update depending on the selected key.
  def handle_key_change(change):
    if change['new'] == custom_annotation_key:
      keyvalue.value = None
    else:
      keyvalue.value = str(df_sample[key.value].iloc[0])

  key.observe(handle_key_change, names='value')

  # Allow the user to leave a text comment as the main value of the annotation.
  comment = widgets.Textarea(
      value='',
      placeholder='Type your comment here',
      description='Comment:',
      disabled=False,
      layout=widgets.Layout(width='80%', height='50px'),
      style={'description_width': 'initial'},
  )

  # Configure the submission button.
  submit_button = widgets.Button(description='Submit annotation', button_style='success')
  output = widgets.Output()

  def on_button_clicked(b):
    params = _format_annotation(sample_id=sample_id, key=key.value, keyvalue=keyvalue.value, comment=comment.value)
    try:
      success = annotation_storage.submit_annotation(
          sample_id=params['sample_id'],
          annotator=params['annotator'],
          key=params['key'],
          value_numeric=params['value_numeric'],
          value_string=params['value_string'],
          comment=params['comment'],
      )
    except Exception as e:
      display(
          HTML(f'''<div class="alert alert-block alert-danger">
                   <b>Warning:</b> Unable to store annotation.
                   <hr><p><pre>{e}</pre></p>
                   </div>'''),
      )
      return()
    with output:
      if success:  # Show the information that was submitted.
        display(
            HTML(f'''<div class="alert alert-block alert-info">
                     Submission successful\n[{annotation_storage.describe()}]</div>'''),
        )
        display(annotation_storage.view_recent_submissions(1))
      else:
        display(
            HTML('''<div class="alert alert-block alert-warning">
                    Annotation not submitted. Please try again.</div>'''),
        )

  submit_button.on_click(on_button_clicked)

  # Display all the widgets.
  display(sample, box1, comment, submit_button, output)


def _format_annotation(sample_id, key, keyvalue, comment):
  """Helper method to clean and reshape info from the widgets and the environment into a dictionary representing the annotation."""
  # Programmatically get the identity of the person running this Terra notebook.
  current_user = os.getenv('OWNER_EMAIL')
  # Also support other environments such as AI Platform Notebooks.
  if current_user is None:
    current_user = socket.gethostname()  # By convention, we prefix the hostname with our username.

  # Check whether the value is string or numeric.
  if keyvalue is None:
    value_numeric = None
    value_string = None
  else:
    try:
      value_numeric = float(keyvalue)  # this will fail if the value is text
      value_string = None
    except ValueError:
      value_numeric = None
      value_string = keyvalue

  # Format into a dictionary.
  params = {
      'sample_id': str(sample_id),
      'annotator': current_user,
      'key': key,
      'value_numeric': value_numeric,
      'value_string': value_string,
      'comment': comment,
  }

  return params
