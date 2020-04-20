"""Methods for collecting and submitting annotations within notebooks."""
import os
import socket
import ipywidgets as widgets
from IPython.display import display
from google.cloud import bigquery


def get_df_sample(sample_info, sample_id):
    """ return a dataframe containing only the row for the indicated sample_id
    """

    df_sample = sample_info[sample_info['sample_id'] == str(sample_id)]
    if 0 == df_sample.shape[0]: df_sample = sample_info.query('sample_id == ' + str(sample_id))
    return df_sample


def display_annotation_collector(sample_info, sample_id):
    """Method to create a gui (set of widgets) through which the user can create an annotation
  Args:
    sample_info: dataframe containing all the samples and data
    sample_id: The selected sample for which the values will be displayed.
  """

    df_sample = get_df_sample(sample_info, sample_id)

    # show the sample ID for this annotation
    sample = widgets.HTML(value=f"For sample <b>{sample_id}</b>")

    # allow the user to pick a key about which to comment
    key = widgets.Dropdown(
        options=sample_info.keys(),
        description='Key:',
        disabled=False,
    )

    # return the sample's value for that key
    valuelabel = widgets.Label(value='Value: ')
    keyvalue = widgets.Label(value=str(df_sample[key.value].iloc[0]))

    box1 = widgets.HBox(
        [key, valuelabel, keyvalue],
        layout=widgets.Layout(width='50%'),
    )

    # have keyvalue auto update depending on the selected key
    def handle_key_change(change):
        keyvalue.value = str(df_sample[key.value].iloc[0])

    key.observe(handle_key_change, names='value')

    # allow the user to leave a text comment
    comment = widgets.Textarea(
        value='',
        placeholder='Type your comment here',
        description=f'Comment:',
        disabled=False,
        layout=widgets.Layout(width='80%', height='50px'),
        style={'description_width': 'initial'},
    )

    # configure submission button
    submit_button = widgets.Button(description="Submit annotation")
    output = widgets.Output()

    def get_key(): return key
    def get_keyvalue(): return keyvalue
    def get_comment(): return comment

    def on_button_clicked(b):
        params = format_annotation(sample_id, [get_key(), get_keyvalue(), get_comment()])
        success = bq_submission(params) # returns boolean True if submission succeeded
        with output:
            if success: # show the information that was submitted
                print('Submission successful:')
                display(view_submissions(1))
            else:
                print('Annotation not submitted. Please try again.\n') # TODO give more information on failure
    submit_button.on_click(on_button_clicked)

    # display everything
    display(sample, box1, comment, submit_button, output)


def format_annotation(sample_id, annotation_data):
    # pull out values from output
    key = annotation_data[0].value
    keyvalue = annotation_data[1].value
    comment = annotation_data[2].value

    # Programmatically get the identity of the person running this Terra notebook.
    USER = os.getenv('OWNER_EMAIL')
    # Also support other environments such as AI Platform Notebooks.
    if USER is None:
        USER = socket.gethostname() # By convention, we prefix the hostname with our username.

    # check whether the value is string or numeric
    if keyvalue is None:
        value_numeric = None
        value_string = None
    else:
        try:
            value_numeric = float(keyvalue)  # this will fail if the value is text
            value_string = None
        except:
            value_numeric = None
            value_string = keyvalue

    # format into a dictionary
    params = {
        'sample_id': str(sample_id),
        'annotator': USER,
        'key': key,
        'value_numeric': value_numeric,
        'value_string': value_string,
        'comment': comment,
    }

    return params


def bq_submission(params, table='uk-biobank-sek-data.ml_results.annotations'):
    """ call a bigquery insert statement to add a row containing annotation information containing
    params (a dict created/formatted by format_annotation)
    """
    # set up biquery client
    bqclient = bigquery.Client(credentials=bigquery.magics.context.credentials)

    # format the insert string
    query_string = '''
INSERT INTO `{table}`
(sample_id, annotator, annotation_timestamp, key, value_numeric, value_string, comment)
VALUES
('{sample_id}', '{annotator}', CURRENT_TIMESTAMP(), '{key}', SAFE_CAST('{value_numeric}' as NUMERIC), '{value_string}', '{comment}')
'''.format(
        table=table,
        sample_id=params['sample_id'],
        annotator=params['annotator'],
        key=params['key'],
        value_numeric=params['value_numeric'],
        value_string=params['value_string'],
        comment=params['comment']
    )

    # submit the insert request
    submission = bqclient.query(query_string)

    # return True if the submission completed TODO test behavior when submission fails
    return submission.done()


def view_submissions(count=10, table='uk-biobank-sek-data.ml_results.annotations'):
    """ view a list of up to [count] most recent submissions from the user
    """

    # set up biquery client
    bqclient = bigquery.Client(credentials=bigquery.magics.context.credentials)

    # format the query string
    query_string = '''SELECT * FROM `{table}`
ORDER BY annotation_timestamp DESC
LIMIT {count}'''.format(
        table=table,
        count=str(count)
    )

    # submit the query and store the result as a dataframe
    df = bqclient.query(query_string).result().to_dataframe()

    return df
