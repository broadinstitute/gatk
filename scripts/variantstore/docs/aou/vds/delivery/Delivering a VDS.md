# Delivering a VDS

## Initial VDS Creation

The Hail VDSes produced by GVS workflows are created within the callset workspace bucket. These are the usual
Terra workspace buckets with names like `gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375`. Once a VDS has been
created, the Variants team works in conjunction with Lee to QA the VDS before we consider delivering it to AoU.

## VDS Delivery to AoU

### Determine the VDS destination path

Once the VDS has passed QA the Variants team can begin the process of delivery to AoU. The designated delivery bucket
for AoU callsets is `gs://prod-drc-broad/`. Within the Variants team we refer to callsets by
their [NATO phonetic alphabet](https://en.wikipedia.org/wiki/NATO_phonetic_alphabet) names (e.g. Charlie, Delta, Echo
etc.), but AoU does not know about these names.

Ask Lee what the delivery path of the VDS should be within the AoU delivery bucket. e.g. for the Delta callset the final
VDS had a source path
of `gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/delta_ai_an_filtered_2022_11_15.vds` and a delivery path
of `gs://prod-drc-broad/v7/wgs/without_aian_prod/vds/aou_srwgs_short_variants_v7_prod.vds`. The source path was named to
make sense to members of the Variants team only, while the delivery path was completely specified by Lee apart from the
standard delivery bucket name.

### Collect callset workspace project info

In order to configure permissions for the VDS transfer we will need both the Google project id and project number. The
project ID can be found in the workspace dashboard at the far left (circled in red in this
example):

![Project number](Callset%20Workspace%20Dashboard.png)

Using our PMI ops account and the example `terra-b3b4392d` project id we can get the project number:

```shell
$ gcloud projects describe terra-b3b4392d --format json | jq -r .projectNumber
496370717638
```

### Open PMI DRC JIRA ticket to enable STS

At this point we will need to ask the PMI DRC to add the appropriate permissions for the STS service account for our
workspace project. This request should include our callset workspace project ID and project number we collected in the
previous step.
For reference [this](https://precisionmedicineinitiative.atlassian.net/browse/PD-8286) was the ticket
we used for this purpose for the Delta callset. For Echo and other future callsets it would be necessary to do at least
the second bullet point in this ticket add the `Storage Object Power User` permission for
the `project-<project number>@storage-transfer-service.iam.gserviceaccount.com` STS service account.

### Grant STS service account permissions on source bucket

**Note:** The following steps required the use of a `@firecloud.org` account and were all performed within the Google
Cloud Console.

For Delta I (Miguel) first went to
the [IAM page for this project](https://console.cloud.google.com/iam-admin/iam?project=terra-b3b4392d)
and granted my  `mcovarr@firecloud.org` principal a `Storage Admin` role in the `terra-b3b4392d`
project:

![Firecloud Storage Admin](Firecloud%20Storage%20Admin.png)

Once my principal had this role I could navigate to `Cloud Storage` and see the workspace bucket:

![Workspace bucket](./Workspace%20Bucket.png)

Clicking into the bucket, and then onto `Permissions`, I assigned the `Storage Object Admin` role to the STS service
account for the callset workspace project:

![Service Account Permissions](./Service%20Account%20Permissions.png)

### Initiate STS transfer

**Note:** The following steps required the use of a `@firecloud.org` account and were all performed within the Google
Cloud Console.

insert notes and George's screenshots here
