import json
import requests
import argparse

from terra_notebook_utils import gs

def get_workspace_name(workspace_id, workspace_name_output, workspace_namespace_output):
    with open(workspace_name_output, "w") as name_output, open(workspace_namespace_output, "w") as namespace_output:
        token = gs.get_access_token()
        # grab the workspace information from rawls
        rawls = 'https://rawls.dsde-prod.broadinstitute.org/api/workspaces/id/{}?fields=workspace.name,workspace.namespace,workspace.googleProject'.format(workspace_id)
        head = {'Authorization': 'Bearer {}'.format(token)}
        response = requests.get(rawls, headers=head)
        ## TODO add an error msg if we get a 400 etc
        response_dict = json.loads(response.text)
        # then extract the googleProject info
        # google_project_id=response_dict['workspace']['googleProject']
        workspace_name=response_dict['workspace']['name']
        name_output.write(f'{workspace_name}\n')
        workspace_namespace=response_dict['workspace']['namespace']
        namespace_output.write(f'{workspace_namespace}\n')
        return (workspace_namespace, workspace_name)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Get workspace information')

    parser.add_argument('--workspace_id', type=str,
                        help='The ID of your workspace that holds your sample data',
                        required=True)

    parser.add_argument('--workspace_name_output', type=str,
                        help='The location to write the workspace name to',
                        required=False)

    parser.add_argument('--workspace_namespace_output', type=str,
                        help='The location to write the workspace namespace to',
                        required=False)


    args = parser.parse_args()


    (workspace_namespace, workspace_name) = get_workspace_name(args.workspace_id, args.workspace_name_output, args.workspace_namespace_output)
