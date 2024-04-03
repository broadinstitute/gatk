# How to Monitor and Request Changes to GCP Quotas

## General Requirements

- an active `@firecloud.org` account
- the name of the workspace Google project where the workflow that you are monitoring is running (e.g. `terra-abcdefgh`)
- the "Service Usage Viewer" Role in the workspace Google project (you can add this to your `@firecloud.org` account using your `@firecloud.org` account at https://console.cloud.google.com/iam-admin/iam?project=[workspace-project-id])

## Monitoring

1. navigate to https://console.cloud.google.com/iam-admin/quotas?project=[workspace-project-id]
1. click on the heading "Current usage percentage" in the table of quotas so that it is sorting the rows in descending usage, which will bring the quotas that are closest to being hit to the top
1. if you are curious about the usage for particular quota in the past, click on the "show usage chart" icon in the far right of that row

## Requesting a Change

1. you will need to upgrade your `@firecloud.org` role to "Service Usage Admin" in the workspace Google project at https://console.cloud.google.com/iam-admin/iam?project=[workspace-project-id] using your `@firecloud.org` account
1. confirm that your `@firecloud.org` has gmail enabled
1. navigate to https://console.cloud.google.com/iam-admin/quotas?project=[workspace-project-id]
1. click on the checkbox at the far left of the row for the quota in question and the "EDIT QUOTAS" button at the top right of the page
1. enter a new value (when in doubt, double the current value) and a request description (i.e. "VCF extraction for newest All Of Us callset") and then click the "NEXT" button
1. fill out your name and click on the "SUBMIT REQUEST" button 
1. you can monitor the status of your request in the "INCREASE REQUESTS" tab at https://console.cloud.google.com/iam-admin/quotas?project=[workspace-project-id]

## Notes

- quota increase requests usually fall into one of two categories
  1. will be automatically granted within minutes of request, examples:
     - Compute Engine API: Persistent Disk Standard (GB)
     - Compute Engine API: VM instances
     - Compute Engine API: CPUs
  1. will need human intervention, which means that, after you make the request, (or even beforehand, if possible) contact our Google Support Account people and ask them to follow up, examples:
     - Compute Engine API: In-use IP addresses
     - BigQuery: Create WriteStreamRequests quota for US regions per minute
- a `@firecloud.org` will not be able to look at the quotas for resources that are in VUMC-managed Google projects (e.g. `aou-genomics-curation-prod`), only Terra-related projects
