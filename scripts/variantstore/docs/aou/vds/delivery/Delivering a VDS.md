# Delivering a VDS

## VDS Creation and QA

The Hail VDSes produced by GVS workflows are created within the callset workspace bucket. These are the usual
Terra workspace buckets with names like `gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375`. Once a VDS has been
created the Variants team will also generate callset statistics using `GvsCallsetStatistics.wdl`. The Variants team
then forwards both the workspace path to the VDS and the output callset statistics TSV to Lee to quality check the VDS.

## VDS Delivery to AoU

### Have Lee determine the VDS delivery path

Once Lee has completed his QA procedure, the Variants team can begin the process of delivering the VDS to AoU. The
designated delivery bucket for AoU callsets is `gs://prod-drc-broad/`. Within the Variants team we refer to callsets
by their [NATO phonetic alphabet](https://en.wikipedia.org/wiki/NATO_phonetic_alphabet) names (e.g. Charlie, Delta, Echo
etc.), but AoU does not use these names.

Lee should determine what the delivery path of the VDS will be within the AoU delivery bucket. e.g. for
the Delta callset the final VDS had a workspace path
of `gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/delta_ai_an_filtered_2022_11_15.vds` and a delivery path
of `gs://prod-drc-broad/v7/wgs/without_aian_prod/vds/aou_srwgs_short_variants_v7_prod.vds`. The source path was named to
make sense to members of the Variants team, while the delivery path was completely specified by Lee, apart from the
standard delivery bucket name.

### Collect callset workspace project info

Both the workspace Google project id and project number are required to configure permissions for VDS delivery. The
project ID can be found in the workspace dashboard at the far left (circled in red):

![Project number](Callset%20Workspace%20Dashboard.png)

Using a PMI ops account and the example `terra-b3b4392d` project id we can get the project number:

```shell
$ gcloud projects describe terra-b3b4392d --format json | jq -r .projectNumber
496370717638
```

### Open PMI DRC JIRA ticket to enable STS

At this point the Variants team needs to ask the PMI DRC to add the appropriate permissions for the STS service account
for our workspace project. This request should include our callset workspace project ID and project number we collected
in the previous step.
For reference [this](https://precisionmedicineinitiative.atlassian.net/browse/PD-8286) was the ticket used for this
purpose for the Delta callset. For Echo and other future callsets it would be necessary to do at least the second bullet
point in this ticket add the `Storage Object Power User` permission for
the `project-<project number>@storage-transfer-service.iam.gserviceaccount.com` STS service account.

### Grant STS service account permissions on source bucket

**Note:** The following steps required the use of a `@firecloud.org` account and were all performed within the Google
Cloud Console.

For the Delta callset, I (Miguel) first went to
the [IAM page for this project](https://console.cloud.google.com/iam-admin/iam?project=terra-b3b4392d)
and granted my  `mcovarr@firecloud.org` principal a `Storage Admin` role in the `terra-b3b4392d` project:

![Firecloud Storage Admin](Firecloud%20Storage%20Admin.png)

Once my principal had this role I could navigate to `Cloud Storage` and see the workspace bucket:

![Workspace bucket](./Workspace%20Bucket.png)

Clicking into the bucket, and then onto `Permissions`, I assigned the `Storage Object Admin` role to the STS service
account for the callset workspace project:

![Service Account Permissions](./Service%20Account%20Permissions.png)

### Perform STS transfer

**Note:** The following steps required the use of a `@firecloud.org` account and were all performed within the Google
Cloud Console.

**Note:** For the Delta callset the Variants
team [delivered two versions of the VDS](https://broadworkbench.atlassian.net/browse/VS-716), one with AI/AN samples and
one without. Since the text and screenshots in this document draw from both deliveries the paths and principals involved
are not consistent. Miguel's principal was used for the version of the VDS without AI/AN and George's principal was used
for the version of the VDS with AI/AN.

#### Go to the Storage Transfer Service (STS) Transfer Jobs Page

From the Google Cloud Console, enter 'Storage Transfer Service' in the search bar and select the 'Data Transfer' result
(circle in red below):

![STS Transfer Jobs Page](./STS%200.png)

This Transfer Jobs page allows for the management of STS jobs for the project. At the top, click the 'Create Transfer
Job' "button":

![STS Create Transfer Job](./STS%201.png)

Clicking this button initiates a 5-step workflow. The first step is to choose the transfer job type, which for our use
case is Google Cloud Storage to Google Cloud Storage:

![STS Get Started](./STS%202.png)

Next choose the source path. This is a GCS URL that can be either a bucket or a directory but the `gs://` prefix must
not be included:

![STS Source](./STS%203.png)

Now choose the destination path. Like the preceding step, this is a GCS URL that can be either a bucket or a directory but the `gs://` prefix must
not be included:

![STS Destination](./STS%204.png)

Take the defaults of "Run once" and "Starting now":

![STS How and when to run](./STS%205.png)

Finally, provide a description to identify the job and choose some conservative options for transfer:

* When to overwrite: Never
* When to delete: Never

![STS Identify job](./STS%206.png)

