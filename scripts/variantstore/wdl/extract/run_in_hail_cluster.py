import argparse
import os
import uuid
from google.cloud import dataproc_v1 as dataproc
from logging import info


def configure_logging():
    import logging
    import sys
    # https://stackoverflow.com/a/14058475
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    root.addHandler(handler)


def wrap(string):
    import re
    return re.sub("\\s{2,}", " ", string).strip()


def run_in_cluster(cluster_name, prefix, contig, account, num_workers, worker_machine_type, region, gcs_project,
                   script_path, vds_url, bed_url, vcf_header_url):

    try:
        cluster_start_cmd = wrap(f"""
        
        hailctl dataproc start 
         --num-workers {num_workers}
         --worker-machine-type {worker_machine_type}
         --region {region}
         --project {gcs_project}
         --service-account {account}
         --num-master-local-ssds 1
         --num-worker-local-ssds 1 
         --max-idle=60m
         --max-age=1440m
         --subnet=projects/{gcs_project}/regions/{region}/subnetworks/subnetwork
         {cluster_name}
         
        """)

        info(f"Starting cluster '{cluster_name}'...")
        info(cluster_start_cmd)
        info(os.popen(cluster_start_cmd).read())

        cluster_client = dataproc.ClusterControllerClient(
            client_options={"api_endpoint": f"{region}-dataproc.googleapis.com:443"}
        )

        for cluster in cluster_client.list_clusters(request={"project_id": gcs_project, "region": region}):
            if cluster.cluster_name == cluster_name:
                cluster_staging_bucket = cluster.config.temp_bucket

                # THIS IS WHERE YOU CALL YOUR SCRIPT AND COPY THE OUTPUT LOCALLY (to get it back into WDL-space)
                submit_cmd = wrap(f"""
                
                gcloud dataproc jobs submit pyspark {script_path}
                 --cluster={cluster_name}
                 --project {gcs_project}
                 --region={region}
                 --account {account}
                 --driver-log-levels root=WARN
                 --
                 --vds_url {vds_url}
                 --bed_url {bed_url}
                 --vcf_header_url {vcf_header_url}
                 --contig {contig}
                 --output_gs_url gs://{cluster_staging_bucket}/{cluster_name}/{prefix}.{contig}.vcf.bgz
                 
                """)

                info("Running: " + submit_cmd)
                os.popen(submit_cmd).read()
                info("Copying results out of staging bucket...")
                staging_cmd = wrap(f"""
                
                gsutil cp -r gs://{cluster_staging_bucket}/{cluster_name}/{prefix}.{contig}.vcf.bgz '{prefix}.{contig}.vcf.bgz'
                
                """)

                info(staging_cmd)
                os.popen(staging_cmd).read()
                ###########
                break
    except Exception as e:
        info(e)
        raise
    finally:
        info(f'Stopping cluster: {cluster_name}')
        delete_cmd = wrap(f"""
            
            gcloud dataproc clusters delete --project {gcs_project} --region {region} --account {account} {cluster_name}
            
        """)

        os.popen(delete_cmd).read()


if __name__ == "__main__":
    configure_logging()

    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Get workspace information')

    parser.add_argument('--cluster-name', type=str, required=True, help='Name of the Hail cluster')
    parser.add_argument('--prefix', type=str, required=True, help='Prefix for output VCF name')
    parser.add_argument('--contig', type=str, required=True, help='Contig to extract')
    parser.add_argument('--account', type=str, required=True, help='GCP account name')
    parser.add_argument('--num-workers', type=str, required=True, help='Number of workers in Hail cluster')
    parser.add_argument('--worker-machine-type', type=str, required=False, help='Dataproc cluster worker machine type')
    parser.add_argument('--region', type=str, required=True, help='GCS region')
    parser.add_argument('--gcs-project', type=str, required=True, help='GCS project')
    parser.add_argument('--script-path', type=str, required=True, help='Path to script to run in Hail cluster')
    parser.add_argument('--vds-url', type=str, required=True, help='VDS URL')
    parser.add_argument('--bed-url', type=str, required=True, help='Bed URL')
    parser.add_argument('--vcf-header-url', type=str, required=True, help='VCF Header URL')


    args = parser.parse_args()

    run_in_cluster(cluster_name=args.cluster_name,
                   prefix=args.prefix,
                   contig=args.contig,
                   account=args.account,
                   num_workers=args.num_workers,
                   worker_machine_type=args.worker_machine_type if args.worker_machine_type else "n1-standard-8",
                   region=args.region,
                   gcs_project=args.gcs_project,
                   script_path=args.script_path,
                   vds_url=args.vds_url,
                   bed_url=args.bed_url,
                   vcf_header_url=args.vcf_header_url
                   )
