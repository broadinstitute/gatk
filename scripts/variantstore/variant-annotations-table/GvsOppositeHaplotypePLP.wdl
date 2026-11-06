version 1.0

import "../wdl/GvsUtils.wdl" as Utils

# For every starting variant in a P/LP-carrying protein-coding gene, find the samples that carry a pathogenic or
# likely pathogenic variant on the *other* haplotype of the same gene (compound heterozygotes, plus the hom-alt cases
# that are trivially on both haplotypes).
#
# Stage A (PrepareShards) reads the BigQuery tables built during scoping and turns them into
#   - a global list of P/LP vids for the genes in scope, with the VCF regions that hold them, and
#   - shards of starting vids, each with the merged VCF regions that hold them.
# PhasePLPSites streams the P/LP regions out of the phased VCF once, over all samples, and records which haplotype
# each carrier has each P/LP variant on. That fixes the carriers, so each ScanShard task then streams only its own
# regions, subset to those carriers, and emits (starting vid, gene, sample) triples.
#
# Only the regions named in the region files are read from the VCF: bcftools -R uses the tabix/CSI index to seek, and
# htslib reads gs:// directly, so the VCF is never localized. Note that sample subsetting does not reduce read cost --
# every record still carries every sample's genotype and is decompressed and tokenized in full.
#
# Keys: VAT vids are minimized and left-aligned, so every VCF record is split (`norm -m-any`, which preserves phase:
# `1|2` becomes `1|0` and `0|1`) and left-aligned (`norm -f`) before it is compared. Indel regions are padded by
# `indel_padding_bp` on both sides because a non-normalized record can sit some distance from its left-aligned
# position. The per-shard stats report how many starting vids matched a VCF record, SNVs and indels separately, which
# is the Q17 join rate in the design doc.
workflow GvsOppositeHaplotypePLP {
    input {
        # The phased joint VCF and its index. Never localized; read remotely by region.
        File phased_vcf
        File? phased_vcf_index

        String project_id
        String dataset
        String partner_vids_table_name = "preupdate_p_lp_gene_partner_vids"
        String plp_vids_table_name = "preupdate_p_lp_vat_entries"

        # Contig to process, named as in the VCF. Intervals on other contigs are ignored.
        String contig = "chr20"
        # Optional BED (0-based half-open) or Picard .interval_list (1-based inclusive) narrowing the *starting*
        # variants. The P/LP variants of every gene touched are always included, wherever they lie in the gene, so
        # the answer for each starting variant inside the intervals is complete.
        File? intervals

        Int vids_per_shard = 50000
        Int indel_padding_bp = 1000
        Int region_merge_gap_bp = 1000

        # Set if the bucket holding the phased VCF is requester pays.
        String? requester_pays_project

        Int plp_disk_gb = 200
        Int shard_disk_gb = 100
        Int shard_cpu = 4
        Int shard_preemptible = 1

        String reference_name = "hg38"
    }

    call Utils.GetToolVersions

    call Utils.GetReference {
        input:
            reference_name = reference_name,
            basic_docker = GetToolVersions.basic_docker,
    }

    call PrepareShards {
        input:
            project_id = project_id,
            dataset = dataset,
            partner_vids_table_name = partner_vids_table_name,
            plp_vids_table_name = plp_vids_table_name,
            contig = contig,
            intervals = intervals,
            vids_per_shard = vids_per_shard,
            indel_padding_bp = indel_padding_bp,
            region_merge_gap_bp = region_merge_gap_bp,
            variants_docker = GetToolVersions.variants_docker,
    }

    call PhasePLPSites {
        input:
            phased_vcf = phased_vcf,
            phased_vcf_index = phased_vcf_index,
            plp_vids = PrepareShards.plp_vids,
            plp_regions = PrepareShards.plp_regions,
            reference = GetReference.reference,
            requester_pays_project = requester_pays_project,
            disk_gb = plp_disk_gb,
            variants_docker = GetToolVersions.variants_docker,
    }

    scatter (shard in PrepareShards.shards) {
        call ScanShard {
            input:
                shard = shard,
                phased_vcf = phased_vcf,
                phased_vcf_index = phased_vcf_index,
                plp_vids = PrepareShards.plp_vids,
                plp_haps = PhasePLPSites.plp_haps,
                reference = GetReference.reference,
                requester_pays_project = requester_pays_project,
                indel_padding_bp = indel_padding_bp,
                cpu = shard_cpu,
                disk_gb = shard_disk_gb,
                preemptible = shard_preemptible,
                variants_docker = GetToolVersions.variants_docker,
        }
    }

    call GatherResults {
        input:
            pairs = ScanShard.pairs,
            undetermined = ScanShard.undetermined,
            unmatched_start_vids = ScanShard.unmatched_start_vids,
            classified_unmatched = ScanShard.classified_unmatched,
            shard_stats = ScanShard.stats,
            prepare_stats = PrepareShards.stats,
            plp_stats = PhasePLPSites.stats,
            variants_docker = GetToolVersions.variants_docker,
    }

    output {
        File pairs = GatherResults.all_pairs
        File undetermined = GatherResults.all_undetermined
        File unmatched_start_vids = GatherResults.all_unmatched_start_vids
        File classified_unmatched = GatherResults.all_classified_unmatched
        File unmatched_plp_vids = PhasePLPSites.unmatched_plp_vids
        File plp_haps = PhasePLPSites.plp_haps
        File summary = GatherResults.summary
    }
}

task PrepareShards {
    input {
        String project_id
        String dataset
        String partner_vids_table_name
        String plp_vids_table_name
        String contig
        File? intervals
        Int vids_per_shard
        Int indel_padding_bp
        Int region_merge_gap_bp
        String variants_docker
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        mkdir shards

        python3 - "~{project_id}" "~{dataset}" "~{partner_vids_table_name}" "~{plp_vids_table_name}" \
            "~{contig}" "~{default='' intervals}" ~{vids_per_shard} ~{indel_padding_bp} ~{region_merge_gap_bp} \
            <<'PYTHON_HEREDOC'
        import json
        import os
        import sys
        import tarfile

        from google.cloud import bigquery

        (project_id, dataset, partner_table, plp_table, contig, intervals_path,
         vids_per_shard, indel_padding, merge_gap) = sys.argv[1:]
        vids_per_shard = int(vids_per_shard)
        indel_padding = int(indel_padding)
        merge_gap = int(merge_gap)


        def read_intervals(path):
            """0-based half-open (start, end) pairs on `contig`, from a BED or a Picard interval list."""
            if not path:
                return [(0, 1 << 40)], 0
            is_interval_list = path.endswith(".interval_list")
            kept, other_contigs = [], 0
            with open(path) as f:
                for line in f:
                    if line.startswith("@"):
                        is_interval_list = True
                        continue
                    if not line.strip() or line.startswith(("#", "track", "browser")):
                        continue
                    fields = line.split()
                    start, end = int(fields[1]), int(fields[2])
                    if is_interval_list:
                        start -= 1
                    if fields[0] != contig:
                        other_contigs += 1
                        continue
                    kept.append((start, end))
            return kept, other_contigs


        def vid_span(vid):
            """1-based inclusive span of the VCF records that could normalize to this vid.

            SNVs and MNPs are not moved by left-alignment, and an unminimized record that contains one still
            overlaps it, which `bcftools view -R` catches by default. Indels can move, so pad them."""
            _, pos, ref, alt = vid.split("-")
            pos = int(pos)
            if len(ref) == len(alt):
                return pos, pos + len(ref) - 1
            return max(1, pos - indel_padding), pos + max(len(ref), len(alt)) + indel_padding


        def write_regions(path, vids):
            merged = []
            for start, end in sorted(vid_span(v) for v in vids):
                if merged and start <= merged[-1][1] + merge_gap + 1:
                    merged[-1][1] = max(merged[-1][1], end)
                else:
                    merged.append([start, end])
            # 1-based inclusive CHROM/BEG/END; bcftools only treats a regions file as 0-based if it is named .bed.
            with open(path, "w") as f:
                for start, end in merged:
                    f.write(f"{contig}\t{start}\t{end}\n")
            return len(merged), sum(end - start + 1 for start, end in merged)


        shard_stats = []


        def write_shard(rows):
            name = f"shard_{len(shard_stats):05d}"
            os.makedirs(name)
            with open(f"{name}/start_vids.tsv", "w") as f:
                for vid, gene_id in rows:
                    f.write(f"{vid}\t{gene_id}\n")
            n_regions, bp = write_regions(f"{name}/start_regions.tsv", {vid for vid, _ in rows})
            with tarfile.open(f"shards/{name}.tar", "w") as tar:
                for member in ("start_vids.tsv", "start_regions.tsv"):
                    tar.add(f"{name}/{member}", arcname=member)
            shard_stats.append({"shard": name, "vid_gene_rows": len(rows), "regions": n_regions, "region_bp": bp})


        intervals, other_contig_intervals = read_intervals(intervals_path)
        if not intervals:
            sys.exit(f"No intervals on {contig} in {intervals_path}")

        client = bigquery.Client(project=project_id)
        contig_param = bigquery.ScalarQueryParameter("contig", "STRING", contig)
        interval_param = bigquery.ArrayQueryParameter("intervals", "STRUCT", [
            bigquery.StructQueryParameter(
                None,
                bigquery.ScalarQueryParameter("bed_start", "INT64", start),
                bigquery.ScalarQueryParameter("bed_end", "INT64", end))
            for start, end in intervals])

        # Starting vids: every partner vid inside the intervals, in position order so shards are contiguous.
        start_sql = f"""
            SELECT DISTINCT p.vid, p.gene_id, p.position
            FROM `{project_id}.{dataset}.{partner_table}` AS p, UNNEST(@intervals) AS i
            WHERE p.contig = @contig
              AND p.position BETWEEN i.bed_start + 1 AND i.bed_end
            ORDER BY p.position, p.vid, p.gene_id
        """
        job = client.query(start_sql, job_config=bigquery.QueryJobConfig(
            query_parameters=[contig_param, interval_param]))

        genes, rows = set(), []
        n_vids = n_in_shard = 0
        last_vid = last_pos = None
        for row in job.result(page_size=100000):
            if row.vid != last_vid:
                # Never split one position across shards, so a record belongs to exactly one shard.
                if n_in_shard >= vids_per_shard and row.position != last_pos:
                    write_shard(rows)
                    rows, n_in_shard = [], 0
                n_in_shard += 1
                n_vids += 1
            rows.append((row.vid, row.gene_id))
            genes.add(row.gene_id)
            last_vid, last_pos = row.vid, row.position
        if rows:
            write_shard(rows)
        if not n_vids:
            sys.exit(f"No partner vids on {contig} inside the intervals")

        # P/LP vids of every gene touched, wherever in the gene they lie.
        plp_sql = f"""
            SELECT DISTINCT vid, gene_id
            FROM `{project_id}.{dataset}.{plp_table}`
            WHERE contig = @contig AND gene_id IN UNNEST(@genes)
            ORDER BY vid, gene_id
        """
        job = client.query(plp_sql, job_config=bigquery.QueryJobConfig(query_parameters=[
            contig_param, bigquery.ArrayQueryParameter("genes", "STRING", sorted(genes))]))
        plp_rows = [(row.vid, row.gene_id) for row in job.result()]
        with open("plp_vids.tsv", "w") as f:
            for vid, gene_id in plp_rows:
                f.write(f"{vid}\t{gene_id}\n")
        n_plp_regions, plp_bp = write_regions("plp_regions.tsv", {vid for vid, _ in plp_rows})

        with open("prepare_stats.json", "w") as f:
            json.dump({
                "contig": contig,
                "intervals": len(intervals),
                "intervals_on_other_contigs_ignored": other_contig_intervals,
                "start_vids": n_vids,
                "genes": len(genes),
                "genes_with_plp": len({g for _, g in plp_rows}),
                "plp_vids": len({v for v, _ in plp_rows}),
                "plp_regions": n_plp_regions,
                "plp_region_bp": plp_bp,
                "shards": shard_stats,
            }, f, indent=2)
        PYTHON_HEREDOC

        cat prepare_stats.json
    >>>
    runtime {
        docker: variants_docker
        memory: "8 GiB"
        disks: "local-disk 50 HDD"
        cpu: 2
        maxRetries: 3
    }
    output {
        Array[File] shards = glob("shards/*.tar")
        File plp_vids = "plp_vids.tsv"
        File plp_regions = "plp_regions.tsv"
        File stats = "prepare_stats.json"
    }
}

task PhasePLPSites {
    input {
        File phased_vcf
        File? phased_vcf_index
        File plp_vids
        File plp_regions
        Reference reference
        String? requester_pays_project
        Int disk_gb
        String variants_docker
    }
    parameter_meta {
        phased_vcf: {
            localization_optional: true
        }
        phased_vcf_index: {
            localization_optional: true
        }
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # htslib reads GCS_OAUTH_TOKEN once, when the file is opened, so a stream that outlives the token fails on its
        # next reconnect. HTS_AUTH_LOCATION is instead re-read whenever the token it last read is within 60 s of its
        # stated expiry. gcloud hands out tokens with at least a few minutes left, so refresh the file every minute
        # and claim a 2 minute lifetime: no token is then used for longer than about 2 minutes after it was minted.
        unset GCS_OAUTH_TOKEN
        export HTS_AUTH_LOCATION="$(mktemp -d)/token.json"
        write_token() {
            set +o xtrace
            printf '{"access_token": "%s", "token_type": "Bearer", "expires_in": 120}\n' \
                "$(gcloud auth application-default print-access-token)" > "${HTS_AUTH_LOCATION}.tmp"
            mv "${HTS_AUTH_LOCATION}.tmp" "${HTS_AUTH_LOCATION}"
            set -o xtrace
        }
        write_token
        ( while true; do sleep 60; write_token; done ) &
        trap "kill $! 2> /dev/null || true" EXIT
        if [[ -n "~{default='' requester_pays_project}" ]]
        then
            export GCS_REQUESTER_PAYS_PROJECT="~{default='' requester_pays_project}"
        fi

        vcf="~{phased_vcf}"
        if [[ -n "~{default='' phased_vcf_index}" ]]
        then
            vcf="${vcf}##idx##~{default='' phased_vcf_index}"
        fi

        bcftools head "${vcf}" > header.txt
        grep '^#CHROM' header.txt | cut -f 10- | tr '\t' '\n' | wc -l > vcf_sample_count.txt

        # All samples: this is the pass that discovers who the carriers are.
        # GT is the only FORMAT field, but INFO carries Number=A fields (AC, AF) that norm would otherwise split on
        # every record, so drop it first.
        bcftools view --threads 4 --regions-file ~{plp_regions} --output-type u "${vcf}" \
            | bcftools annotate --remove INFO --output-type u \
            | bcftools norm --multiallelics -any --output-type u 2> split.log \
            | bcftools norm --fasta-ref ~{reference.reference_fasta} --check-ref w \
                --output-type b --output plp.bcf 2> realign.log
        cat split.log realign.log

        python3 - ~{plp_vids} plp.bcf vcf_sample_count.txt split.log realign.log <<'PYTHON_HEREDOC'
        import gzip
        import json
        import subprocess
        import sys

        plp_vids_path, bcf, sample_count_path, split_log, realign_log = sys.argv[1:]


        def norm_counts(path):
            """The 'Lines total/split/joined/realigned/...' summary bcftools norm writes to stderr, as a dict."""
            for line in open(path):
                if line.startswith("Lines"):
                    _, names, values = line.split()
                    return dict(zip(names.rstrip(":").split("/"), map(int, values.split("/"))))
            return {}


        def hap_of(gt):
            """'1', '2', '12' (both), 'U' (unphased het), 'H' (haploid alt) or None (no alt allele)."""
            alleles = gt.replace("|", "/").split("/")
            alt = [a not in ("0", ".") for a in alleles]
            if not any(alt):
                return None
            if len(alleles) == 1:
                return "H"
            if all(alt):
                return "12"
            if "|" not in gt:
                return "U"
            return "1" if alt[0] else "2"


        def query(fmt, *args):
            proc = subprocess.Popen(["bcftools", "query", *args, "-f", fmt, bcf], stdout=subprocess.PIPE, text=True)
            for line in proc.stdout:
                yield line.rstrip("\n").split("\t")
            if proc.wait() != 0:
                sys.exit(f"bcftools query failed with exit code {proc.returncode}")


        def key(chrom, pos, ref, alt):
            return f"{chrom.removeprefix('chr')}-{pos}-{ref}-{alt}"


        plp = {line.split("\t")[0] for line in open(plp_vids_path)}
        found = {key(*fields) for fields in query("%CHROM\\t%POS\\t%REF\\t%ALT\\n")} & plp

        counts = {"1": 0, "2": 0, "12": 0, "U": 0, "H": 0}
        carriers_per_vid = {}
        with gzip.open("plp_haps.tsv.gz", "wt") as out:
            for chrom, pos, ref, alt, sample, gt in query("[%CHROM\\t%POS\\t%REF\\t%ALT\\t%SAMPLE\\t%GT\\n]",
                                                          "-i", 'GT="alt"'):
                vid = key(chrom, pos, ref, alt)
                hap = hap_of(gt) if vid in plp else None
                if hap is None:
                    continue
                counts[hap] += 1
                carriers_per_vid[vid] = carriers_per_vid.get(vid, 0) + 1
                out.write(f"{vid}\t{sample}\t{hap}\n")

        with open("unmatched_plp_vids.tsv", "w") as f:
            for vid in sorted(plp - found):
                f.write(f"{vid}\n")


        def is_indel(vid):
            _, _, ref, alt = vid.split("-")
            return len(ref) != len(alt)


        split, realign = norm_counts(split_log), norm_counts(realign_log)
        stats = {
            "vcf_samples": int(open(sample_count_path).read()),
            # Raw VCF records read, how many were multi-allelic, and how many split records left-alignment or
            # trimming changed. A nonzero `records_realigned` means the VCF was not already normalized.
            "records_read": split.get("total", 0),
            "records_multiallelic": split.get("split", 0),
            "records_realigned": realign.get("realigned", 0),
            "plp_vids": len(plp),
            "plp_snvs": sum(1 for v in plp if not is_indel(v)),
            "plp_indels": sum(1 for v in plp if is_indel(v)),
            "plp_vids_found": len(found),
            "plp_snvs_unmatched": sum(1 for v in plp - found if not is_indel(v)),
            "plp_indels_unmatched": sum(1 for v in plp - found if is_indel(v)),
            "plp_vids_with_carriers": len(carriers_per_vid),
            "max_carriers_of_one_plp_vid": max(carriers_per_vid.values(), default=0),
            "carrier_calls_by_haplotype": counts,
            "ref_mismatches": sum(1 for line in open(realign_log) if line.startswith("REF_MISMATCH")),
        }
        with open("plp_stats.json", "w") as f:
            json.dump(stats, f, indent=2)
        PYTHON_HEREDOC

        cat plp_stats.json
    >>>
    runtime {
        docker: variants_docker
        memory: "8 GiB"
        disks: "local-disk ~{disk_gb} SSD"
        cpu: 4
        maxRetries: 3
    }
    output {
        File plp_haps = "plp_haps.tsv.gz"
        File unmatched_plp_vids = "unmatched_plp_vids.tsv"
        File stats = "plp_stats.json"
    }
}

task ScanShard {
    input {
        File shard
        File phased_vcf
        File? phased_vcf_index
        File plp_vids
        File plp_haps
        Reference reference
        String? requester_pays_project
        Int indel_padding_bp
        Int cpu
        Int disk_gb
        Int preemptible
        String variants_docker
    }
    parameter_meta {
        phased_vcf: {
            localization_optional: true
        }
        phased_vcf_index: {
            localization_optional: true
        }
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # See PhasePLPSites for why the token goes through HTS_AUTH_LOCATION rather than GCS_OAUTH_TOKEN.
        unset GCS_OAUTH_TOKEN
        export HTS_AUTH_LOCATION="$(mktemp -d)/token.json"
        write_token() {
            set +o xtrace
            printf '{"access_token": "%s", "token_type": "Bearer", "expires_in": 120}\n' \
                "$(gcloud auth application-default print-access-token)" > "${HTS_AUTH_LOCATION}.tmp"
            mv "${HTS_AUTH_LOCATION}.tmp" "${HTS_AUTH_LOCATION}"
            set -o xtrace
        }
        write_token
        ( while true; do sleep 60; write_token; done ) &
        trap "kill $! 2> /dev/null || true" EXIT
        if [[ -n "~{default='' requester_pays_project}" ]]
        then
            export GCS_REQUESTER_PAYS_PROJECT="~{default='' requester_pays_project}"
        fi

        vcf="~{phased_vcf}"
        if [[ -n "~{default='' phased_vcf_index}" ]]
        then
            vcf="${vcf}##idx##~{default='' phased_vcf_index}"
        fi

        tar xf ~{shard}

        # Carriers: samples with a P/LP call in any gene this shard touches. Only they can produce an answer here.
        cut -f 2 start_vids.tsv | sort -u > shard_genes.txt
        awk -F '\t' 'NR == FNR { genes[$1]; next } ($2 in genes) { print $1 }' shard_genes.txt ~{plp_vids} \
            | sort -u > shard_plp_vids.txt
        gzip -dc ~{plp_haps} \
            | awk -F '\t' 'NR == FNR { vids[$1]; next } ($1 in vids) { print $2 }' shard_plp_vids.txt - \
            | sort -u > carriers.txt
        wc -l shard_genes.txt shard_plp_vids.txt carriers.txt

        if [[ -s carriers.txt ]]
        then
            bcftools view --threads ~{cpu} --regions-file start_regions.tsv --samples-file carriers.txt --no-update \
                --output-type u "${vcf}" \
                | bcftools annotate --remove INFO --output-type u \
                | bcftools norm --multiallelics -any --output-type u 2> split.log \
                | bcftools norm --fasta-ref ~{reference.reference_fasta} --check-ref w \
                    --output-type b --output start.bcf 2> realign.log
            cat split.log realign.log
        else
            # Nothing to pair against, so skip the VCF entirely. The join-rate stats are then absent for this shard.
            touch split.log realign.log
        fi

        python3 - start_vids.tsv ~{plp_vids} ~{plp_haps} start.bcf split.log realign.log ~{indel_padding_bp} \
            <<'PYTHON_HEREDOC'
        import bisect
        import gzip
        import json
        import os
        import subprocess
        import sys
        from collections import defaultdict

        (start_vids_path, plp_vids_path, plp_haps_path, bcf, split_log, realign_log,
         indel_padding) = sys.argv[1:]
        indel_padding = int(indel_padding)


        def norm_counts(path):
            """The 'Lines total/split/joined/realigned/...' summary bcftools norm writes to stderr, as a dict."""
            for line in open(path):
                if line.startswith("Lines"):
                    _, names, values = line.split()
                    return dict(zip(names.rstrip(":").split("/"), map(int, values.split("/"))))
            return {}


        def hap_of(gt):
            """'1', '2', '12' (both), 'U' (unphased het), 'H' (haploid alt) or None (no alt allele)."""
            alleles = gt.replace("|", "/").split("/")
            alt = [a not in ("0", ".") for a in alleles]
            if not any(alt):
                return None
            if len(alleles) == 1:
                return "H"
            if all(alt):
                return "12"
            if "|" not in gt:
                return "U"
            return "1" if alt[0] else "2"


        def relation(p_hap, v_hap):
            """How a P/LP call and a starting-variant call in one sample relate.

            'phased' when they are on opposite haplotypes, 'homalt' when either is on both (so the other is
            necessarily opposite it, phase or no phase), 'undetermined' when an unphased het leaves it open, else
            None (same haplotype)."""
            if p_hap == "12" or v_hap == "12":
                return "homalt"
            if p_hap == "U" or v_hap == "U":
                return "undetermined"
            return "phased" if p_hap != v_hap else None


        def query(fmt, *args):
            proc = subprocess.Popen(["bcftools", "query", *args, "-f", fmt, bcf], stdout=subprocess.PIPE, text=True)
            for line in proc.stdout:
                yield line.rstrip("\n").split("\t")
            if proc.wait() != 0:
                sys.exit(f"bcftools query failed with exit code {proc.returncode}")


        def key(chrom, pos, ref, alt):
            return f"{chrom.removeprefix('chr')}-{pos}-{ref}-{alt}"


        def is_indel(vid):
            _, _, ref, alt = vid.split("-")
            return len(ref) != len(alt)


        start_genes = defaultdict(list)
        for line in open(start_vids_path):
            vid, gene_id = line.rstrip("\n").split("\t")
            start_genes[vid].append(gene_id)
        shard_genes = {g for genes in start_genes.values() for g in genes}

        plp_genes = defaultdict(list)
        for line in open(plp_vids_path):
            vid, gene_id = line.rstrip("\n").split("\t")
            if gene_id in shard_genes:
                plp_genes[vid].append(gene_id)

        # (gene, sample) -> [(P/LP vid, haplotype)]
        plp_calls = defaultdict(list)
        with gzip.open(plp_haps_path, "rt") as f:
            for line in f:
                vid, sample, hap = line.rstrip("\n").split("\t")
                for gene_id in plp_genes.get(vid, ()):
                    plp_calls[(gene_id, sample)].append((vid, hap))

        stats = defaultdict(int)
        have_bcf = os.path.exists(bcf)
        site_records = list(query("%CHROM\\t%POS\\t%REF\\t%ALT\\n")) if have_bcf else []
        found = {key(*fields) for fields in site_records} & start_genes.keys()
        unmatched = sorted(start_genes.keys() - found) if have_bcf else []
        with open("unmatched_start_vids.tsv", "w") as f:
            for vid in unmatched:
                stats["start_indels_unmatched" if is_indel(vid) else "start_snvs_unmatched"] += 1
                f.write(f"{vid}\n")
        stats["start_vids"] = len(start_genes)
        stats["start_snvs"] = sum(1 for v in start_genes if not is_indel(v))
        stats["start_indels"] = sum(1 for v in start_genes if is_indel(v))
        stats["start_vids_found"] = len(found)
        split, realign = norm_counts(split_log), norm_counts(realign_log)
        stats["records_read"] = split.get("total", 0)
        stats["records_multiallelic"] = split.get("split", 0)
        stats["records_realigned"] = realign.get("realigned", 0)

        seen = set()
        with gzip.open("pairs.tsv.gz", "wt") as pairs, gzip.open("undetermined.tsv.gz", "wt") as undetermined:
            if have_bcf:
                for chrom, pos, ref, alt, sample, gt in query("[%CHROM\\t%POS\\t%REF\\t%ALT\\t%SAMPLE\\t%GT\\n]",
                                                              "-i", 'GT="alt"'):
                    vid = key(chrom, pos, ref, alt)
                    v_hap = hap_of(gt) if vid in start_genes else None
                    if v_hap is None:
                        continue
                    stats["start_calls"] += 1
                    if v_hap == "H":
                        stats["haploid_calls_skipped"] += 1
                        continue
                    for gene_id in start_genes[vid]:
                        if (vid, gene_id, sample) in seen:
                            stats["duplicate_calls_after_normalization"] += 1
                            continue
                        seen.add((vid, gene_id, sample))
                        relations = defaultdict(list)
                        for plp_vid, p_hap in plp_calls.get((gene_id, sample), ()):
                            # A variant is not its own partner, even when hom-alt. See open question 6.
                            if plp_vid != vid:
                                relations[relation(p_hap, v_hap)].append(plp_vid)
                        if relations["phased"] or relations["homalt"]:
                            basis = "phased" if relations["phased"] else "homalt"
                            stats[f"pairs_{basis}"] += 1
                            plps = ",".join(sorted(set(relations["phased"] + relations["homalt"])))
                            pairs.write(f"{vid}\t{gene_id}\t{sample}\t{basis}\t{plps}\n")
                        elif relations["undetermined"]:
                            stats["pairs_undetermined"] += 1
                            plps = ",".join(sorted(set(relations["undetermined"])))
                            undetermined.write(f"{vid}\t{gene_id}\t{sample}\t{plps}\n")

        stats["ref_mismatches"] = sum(1 for line in open(realign_log) if line.startswith("REF_MISMATCH"))
        stats["scanned"] = have_bcf

        # Why each unmatched vid did not match, judged against the shard's own normalized sites. start.bcf keeps
        # every site in the shard's regions even where no carrier has an alt allele, so this needs no further VCF
        # reads.
        by_pos = defaultdict(list)
        indels_by_pos = []
        for _, pos, ref, alt in site_records:
            by_pos[int(pos)].append(f"{ref}>{alt}")
            if len(ref) != len(alt):
                indels_by_pos.append((int(pos), len(alt) - len(ref)))
        indels_by_pos.sort()


        def classify(vid):
            """'other_alleles_at_position' when the VCF has a record at the vid's position with other alleles;
            for indels, 'same_length_indel_nearby' when an indel of the same length change sits within the padding
            (indel length survives normalization, so this suggests a representation difference left-alignment did
            not resolve); otherwise 'no_record_at_position'."""
            _, pos, ref, alt = vid.split("-")
            pos = int(pos)
            if by_pos.get(pos):
                return "other_alleles_at_position", by_pos[pos]
            if len(ref) != len(alt):
                lo = bisect.bisect_left(indels_by_pos, (pos - indel_padding, -(1 << 40)))
                hi = bisect.bisect_right(indels_by_pos, (pos + indel_padding, 1 << 40))
                nearby = [str(p) for p, ilen in indels_by_pos[lo:hi] if ilen == len(alt) - len(ref)]
                if nearby:
                    return "same_length_indel_nearby", nearby
            return "no_record_at_position", []


        with open("classified_unmatched.tsv", "w") as out:
            for vid in unmatched:
                kind = "indel" if is_indel(vid) else "snv"
                category, detail = classify(vid)
                stats[f"unmatched_{kind}_{category}"] += 1
                out.write(f"{vid}\t{kind}\t{category}\t{','.join(detail[:10])}\n")

        with open("shard_stats.json", "w") as f:
            json.dump(stats, f, indent=2)
        PYTHON_HEREDOC

        cat shard_stats.json
    >>>
    runtime {
        docker: variants_docker
        memory: "16 GiB"
        disks: "local-disk ~{disk_gb} SSD"
        cpu: cpu
        preemptible: preemptible
        maxRetries: 3
    }
    output {
        File pairs = "pairs.tsv.gz"
        File undetermined = "undetermined.tsv.gz"
        File unmatched_start_vids = "unmatched_start_vids.tsv"
        File classified_unmatched = "classified_unmatched.tsv"
        File stats = "shard_stats.json"
    }
}

task GatherResults {
    input {
        Array[File] pairs
        Array[File] undetermined
        Array[File] unmatched_start_vids
        Array[File] classified_unmatched
        Array[File] shard_stats
        File prepare_stats
        File plp_stats
        String variants_docker
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # Concatenated gzip members are themselves a valid gzip file.
        printf 'vid\tgene_id\tsample\tbasis\tplp_vids\n' | gzip > pairs.tsv.gz
        cat ~{sep=' ' pairs} >> pairs.tsv.gz
        printf 'vid\tgene_id\tsample\tplp_vids\n' | gzip > undetermined.tsv.gz
        cat ~{sep=' ' undetermined} >> undetermined.tsv.gz
        cat ~{sep=' ' unmatched_start_vids} | gzip > unmatched_start_vids.tsv.gz
        { printf 'vid\ttype\tcategory\tnormalized_records_at_position_or_nearby_indel_positions\n'
          cat ~{sep=' ' classified_unmatched}; } | gzip > classified_unmatched.tsv.gz

        python3 - ~{prepare_stats} ~{plp_stats} ~{sep=' ' shard_stats} <<'PYTHON_HEREDOC'
        import json
        import sys
        from collections import Counter

        prepare = json.load(open(sys.argv[1]))
        plp = json.load(open(sys.argv[2]))
        shards = Counter()
        for path in sys.argv[3:]:
            shards.update(json.load(open(path)))

        prepare.pop("shards")


        def rate(numerator, denominator):
            return numerator / denominator if denominator else None


        summary = {
            "prepare": prepare,
            "plp": plp,
            "shards": dict(sorted(shards.items())),
            "start_vid_match_rate": rate(shards["start_vids_found"], shards["start_vids"]),
            "start_snv_match_rate": rate(shards["start_snvs"] - shards["start_snvs_unmatched"], shards["start_snvs"]),
            "start_indel_match_rate": rate(shards["start_indels"] - shards["start_indels_unmatched"],
                                           shards["start_indels"]),
        }
        with open("summary.json", "w") as f:
            json.dump(summary, f, indent=2)
        PYTHON_HEREDOC

        cat summary.json
    >>>
    runtime {
        docker: variants_docker
        memory: "4 GiB"
        disks: "local-disk 100 HDD"
        cpu: 1
        maxRetries: 3
    }
    output {
        File all_pairs = "pairs.tsv.gz"
        File all_undetermined = "undetermined.tsv.gz"
        File all_unmatched_start_vids = "unmatched_start_vids.tsv.gz"
        File all_classified_unmatched = "classified_unmatched.tsv.gz"
        File summary = "summary.json"
    }
}
