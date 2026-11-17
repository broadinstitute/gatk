version 1.0

import "../wdl/GvsUtils.wdl" as Utils

# For the starting vids that GvsOppositeHaplotypePLP could not find in the phased VCF, estimate from the VAT whether
# the row filters that produced that VCF would have dropped them.
#
# The phased VCF was made from the same VDS by warp's all_of_us/phasing/AoU_VCF_v9.py, which drops a VDS row (an
# unsplit site) if any of these hold:
#   1. it has more than 31 alt alleles
#   2. no alt allele at the site has AC >= 2
#   3. call rate < 0.95
#   4. mean GQ < 30, or mean summed AD over carriers < 14
#   5. the site carries LowQual, NO_HQ_GENOTYPES or ExcessHet
# The VAT supports proxies for 1-3 only. 4 needs per-genotype GQ and AD, which the VAT does not hold. 5 never applies,
# because the VAT already drops filtered sites.
#
# The proxies are approximate, and the direction of the error matters when reading the output:
#   - The VAT computes AC and AN after setting FT-failing genotypes to no-call (hail_create_vat_inputs.py), while the
#     phasing script counts every genotype. So VAT AC and AN are lower bounds, and the proxies for 2 and 3 can flag a
#     site the phasing filter kept, but cannot miss one it would have dropped.
#   - A VDS site is approximated by the VAT vids at the same position. That is exact for SNVs, but left-alignment can
#     move an indel away from its VDS site, so the sibling alleles of an indel, and of anything sharing a position with
#     one, are approximate.
# Because of the first point, the summary also reports the same proxies over the *matched* starting vids, which were
# by definition kept. The rate at which matched vids are flagged measures how far the proxies over-call; a high
# unmatched rate means something only when set against it.
#
# Pass the same contig and intervals as the GvsOppositeHaplotypePLP run, so that the matched set computed here (all
# starting vids minus the unmatched ones) is the one that run used.
workflow GvsExplainUnmatchedPhasingVids {
    input {
        # The classified_unmatched output of GvsOppositeHaplotypePLP.
        File classified_unmatched

        String project_id
        String dataset
        String vat_table_name
        String partner_vids_table_name = "preupdate_p_lp_gene_partner_vids"

        String contig = "chr20"
        File? intervals

        # Samples in the VDS the phasing filter ran on, for the call-rate proxy. If unset, the largest VAT AN among
        # the starting vids, halved, is used instead; that is a lower bound, so call rates come out slightly high.
        Int? vds_sample_count

        # Thresholds from AoU_VCF_v9.py.
        Int max_alt_alleles = 31
        Int min_site_max_ac = 2
        Float min_call_rate = 0.95

        # Positions per detail query. Each query reads the VAT between its first and last position, so on a VAT
        # clustered on (contig, position) the queries together read about one pass over the region. On a VAT clustered
        # on contig alone, every query scans the whole contig.
        Int positions_per_query = 100000
    }

    call Utils.GetToolVersions

    call ExplainUnmatched {
        input:
            classified_unmatched = classified_unmatched,
            project_id = project_id,
            dataset = dataset,
            vat_table_name = vat_table_name,
            partner_vids_table_name = partner_vids_table_name,
            contig = contig,
            intervals = intervals,
            vds_sample_count = vds_sample_count,
            max_alt_alleles = max_alt_alleles,
            min_site_max_ac = min_site_max_ac,
            min_call_rate = min_call_rate,
            positions_per_query = positions_per_query,
            variants_docker = GetToolVersions.variants_docker,
    }

    output {
        File classification = ExplainUnmatched.classification
        File summary = ExplainUnmatched.summary
    }
}

task ExplainUnmatched {
    input {
        File classified_unmatched
        String project_id
        String dataset
        String vat_table_name
        String partner_vids_table_name
        String contig
        File? intervals
        Int? vds_sample_count
        Int max_alt_alleles
        Int min_site_max_ac
        Float min_call_rate
        Int positions_per_query
        String variants_docker
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        python3 - "~{classified_unmatched}" "~{project_id}" "~{dataset}" "~{vat_table_name}" \
            "~{partner_vids_table_name}" "~{contig}" "~{default='' intervals}" "~{default='' vds_sample_count}" \
            ~{max_alt_alleles} ~{min_site_max_ac} ~{min_call_rate} ~{positions_per_query} \
            <<'PYTHON_HEREDOC'
        import gzip
        import json
        import sys
        from collections import Counter

        from google.cloud import bigquery

        (classified_path, project_id, dataset, vat_table, partner_table, contig, intervals_path, sample_count,
         max_alt_alleles, min_site_max_ac, min_call_rate, positions_per_query) = sys.argv[1:]
        sample_count = float(sample_count) if sample_count else None
        max_alt_alleles = int(max_alt_alleles)
        min_site_max_ac = int(min_site_max_ac)
        min_call_rate = float(min_call_rate)
        positions_per_query = int(positions_per_query)


        def read_intervals(path):
            """0-based half-open (start, end) pairs on `contig`, from a BED or a Picard interval list."""
            if not path:
                return [(0, 1 << 40)]
            is_interval_list = path.endswith(".interval_list")
            kept = []
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
                    if fields[0] == contig:
                        kept.append((start, end))
            return kept


        unmatched = {}
        with gzip.open(classified_path, "rt") as f:
            next(f)
            for line in f:
                vid, kind, category = line.rstrip("\n").split("\t")[:3]
                unmatched[vid] = (kind, category)
        chrom = contig[3:] if contig.startswith("chr") else contig
        elsewhere = [vid for vid in unmatched if vid.split("-")[0] != chrom]
        if elsewhere:
            sys.exit(f"{len(elsewhere)} unmatched vids are not on {contig}, e.g. {elsewhere[0]}")

        intervals = read_intervals(intervals_path)
        if not intervals:
            sys.exit(f"No intervals on {contig} in {intervals_path}")

        client = bigquery.Client(project=project_id)
        bytes_processed = 0


        def run(sql, params):
            global bytes_processed
            job = client.query(sql, job_config=bigquery.QueryJobConfig(query_parameters=params))
            rows = list(job.result(page_size=100000))
            bytes_processed += job.total_bytes_processed or 0
            return rows


        def flagged_sql(position_filter):
            """CTEs ending in `flagged`: one row per VAT vid passing `position_filter`, with the site proxies.

            A site is the set of VAT vids at one position, so `position_filter` must keep whole positions."""
            return f"""
                vat_vids AS (
                  SELECT vid, ANY_VALUE(position) AS position,
                         ANY_VALUE(gvs_all_ac) AS ac, ANY_VALUE(gvs_all_an) AS an
                  FROM `{project_id}.{dataset}.{vat_table}`
                  WHERE contig = @contig AND {position_filter}
                  GROUP BY vid
                ),
                sites AS (
                  SELECT position, COUNT(*) AS site_alleles, MAX(ac) AS site_max_ac
                  FROM vat_vids
                  GROUP BY position
                ),
                n AS (
                  SELECT COALESCE(@sample_count, MAX(an) / 2) AS sample_count FROM vat_vids
                ),
                flagged AS (
                  SELECT v.vid, v.ac, v.an, s.site_alleles, s.site_max_ac, n.sample_count,
                         s.site_max_ac < @min_site_max_ac AS fails_max_ac,
                         v.an / (2 * n.sample_count) < @min_call_rate AS fails_call_rate,
                         s.site_alleles > @max_alt_alleles AS fails_alt_max
                  FROM vat_vids AS v
                  JOIN sites AS s USING (position)
                  CROSS JOIN n
                )"""


        def common_params(n_samples):
            return [
                bigquery.ScalarQueryParameter("contig", "STRING", contig),
                bigquery.ScalarQueryParameter("sample_count", "FLOAT64", n_samples),
                bigquery.ScalarQueryParameter("min_site_max_ac", "INT64", min_site_max_ac),
                bigquery.ScalarQueryParameter("min_call_rate", "FLOAT64", min_call_rate),
                bigquery.ScalarQueryParameter("max_alt_alleles", "INT64", max_alt_alleles),
            ]


        flag_names = (f"vat_site_max_ac_lt_{min_site_max_ac}",
                      f"vat_call_rate_lt_{min_call_rate:g}",
                      f"vat_site_alleles_gt_{max_alt_alleles}")
        not_flagged = {"none_vat_checkable", "vat_counts_missing", "not_in_vat"}


        def reason(fails_max_ac, fails_call_rate, fails_alt_max):
            flags = (fails_max_ac, fails_call_rate, fails_alt_max)
            if any(flag is None for flag in flags):
                return "vat_counts_missing"
            return "+".join(name for name, flag in zip(flag_names, flags) if flag) or "none_vat_checkable"


        # Control: the proxies over every starting vid of the original run, aggregated. This also resolves the sample
        # count when it was not supplied, from all starting vids rather than from one detail chunk.
        lo = min(start for start, _ in intervals) + 1
        hi = max(end for _, end in intervals)
        control_sql = f"""
            WITH {flagged_sql("position BETWEEN @lo AND @hi")},
            start_vids AS (
              SELECT DISTINCT p.vid
              FROM `{project_id}.{dataset}.{partner_table}` AS p, UNNEST(@intervals) AS i
              WHERE p.contig = @contig AND p.position BETWEEN i.bed_start + 1 AND i.bed_end
            )
            SELECT LENGTH(SPLIT(vid, '-')[OFFSET(2)]) != LENGTH(SPLIT(vid, '-')[OFFSET(3)]) AS is_indel,
                   fails_max_ac, fails_call_rate, fails_alt_max,
                   COUNT(*) AS vids, ANY_VALUE(sample_count) AS sample_count
            FROM start_vids JOIN flagged USING (vid)
            GROUP BY 1, 2, 3, 4
        """
        control_rows = run(control_sql, common_params(sample_count) + [
            bigquery.ScalarQueryParameter("lo", "INT64", lo),
            bigquery.ScalarQueryParameter("hi", "INT64", hi),
            bigquery.ArrayQueryParameter("intervals", "STRUCT", [
                bigquery.StructQueryParameter(
                    None,
                    bigquery.ScalarQueryParameter("bed_start", "INT64", start),
                    bigquery.ScalarQueryParameter("bed_end", "INT64", end))
                for start, end in intervals]),
        ])
        if not control_rows:
            sys.exit(f"No starting vids on {contig} inside the intervals")
        sample_count_source = "input" if sample_count else "max VAT AN / 2 over the starting vids"
        sample_count = control_rows[0].sample_count
        all_start = Counter()
        for row in control_rows:
            kind = "indel" if row.is_indel else "snv"
            all_start[(kind, reason(row.fails_max_ac, row.fails_call_rate, row.fails_alt_max))] += row.vids

        # Detail: the proxies for each unmatched vid, querying only the positions that hold one.
        detail = {}
        positions = sorted({int(vid.split("-")[1]) for vid in unmatched})
        detail_sql = f"""
            WITH {flagged_sql("position BETWEEN @chunk_lo AND @chunk_hi AND position IN UNNEST(@positions)")}
            SELECT vid, ac, an, site_alleles, site_max_ac, fails_max_ac, fails_call_rate, fails_alt_max
            FROM flagged
        """
        for i in range(0, len(positions), positions_per_query):
            chunk = positions[i:i + positions_per_query]
            # The range is redundant with the list, but it is what lets a position-clustered VAT prune blocks.
            for row in run(detail_sql, common_params(sample_count) + [
                    bigquery.ScalarQueryParameter("chunk_lo", "INT64", chunk[0]),
                    bigquery.ScalarQueryParameter("chunk_hi", "INT64", chunk[-1]),
                    bigquery.ArrayQueryParameter("positions", "INT64", chunk)]):
                if row.vid in unmatched:
                    detail[row.vid] = row

        unmatched_counts = Counter()
        with gzip.open("unmatched_filter_classification.tsv.gz", "wt") as out:
            out.write("vid\ttype\tmatch_category\tvat_ac\tvat_an\tvat_call_rate\tvat_site_alleles\t"
                      "vat_site_max_ac\treason\n")
            for vid, (kind, category) in sorted(unmatched.items()):
                row = detail.get(vid)
                if row is None:
                    fields, why = [""] * 5, "not_in_vat"
                else:
                    call_rate = None if row.an is None else f"{row.an / (2 * sample_count):.4f}"
                    fields = ["" if value is None else str(value)
                              for value in (row.ac, row.an, call_rate, row.site_alleles, row.site_max_ac)]
                    why = reason(row.fails_max_ac, row.fails_call_rate, row.fails_alt_max)
                unmatched_counts[(kind, category, why)] += 1
                out.write("\t".join([vid, kind, category] + fields + [why]) + "\n")

        # Matched = all starting vids minus the unmatched ones, by type and reason. The control only sees vids that are
        # in the VAT, so leave out the not_in_vat ones rather than drive a count negative.
        unmatched_by_reason = Counter()
        for (kind, _, why), n in unmatched_counts.items():
            unmatched_by_reason[(kind, why)] += n
        matched = {key: all_start[key] - unmatched_by_reason[key]
                   for key in all_start.keys() | unmatched_by_reason.keys() if key[1] != "not_in_vat"}
        negative = sorted(f"{kind}/{why}" for (kind, why), n in matched.items() if n < 0)


        def nest(counter, depth):
            out = {}
            for key, n in sorted(counter.items()):
                node = out
                for part in key[:depth - 1]:
                    node = node.setdefault(part, {})
                node[key[depth - 1]] = n
            return out


        def flagged_fraction(counter, kind):
            total = sum(n for (k, _), n in counter.items() if k == kind)
            flagged = sum(n for (k, why), n in counter.items() if k == kind and why not in not_flagged)
            return flagged / total if total else None


        summary = {
            "settings": {
                "contig": contig,
                "intervals": len(intervals),
                "sample_count": sample_count,
                "sample_count_source": sample_count_source,
                "max_alt_alleles": max_alt_alleles,
                "min_site_max_ac": min_site_max_ac,
                "min_call_rate": min_call_rate,
                "bytes_processed": bytes_processed,
            },
            "flagged_fraction": {
                kind: {
                    "unmatched": flagged_fraction(unmatched_by_reason, kind),
                    "matched": flagged_fraction(Counter(matched), kind),
                }
                for kind in ("snv", "indel")
            },
            "unmatched": nest(unmatched_counts, 3),
            "matched": nest(Counter(matched), 2),
            "all_start_vids": nest(all_start, 2),
        }
        if negative:
            summary["warning"] = ("More unmatched than starting vids for " + ", ".join(negative) +
                                  ". The contig or intervals probably differ from the GvsOppositeHaplotypePLP run.")
        with open("unmatched_filter_summary.json", "w") as f:
            json.dump(summary, f, indent=2)
        PYTHON_HEREDOC

        cat unmatched_filter_summary.json
    >>>
    runtime {
        docker: variants_docker
        memory: "8 GiB"
        disks: "local-disk 50 HDD"
        cpu: 2
        maxRetries: 3
    }
    output {
        File classification = "unmatched_filter_classification.tsv.gz"
        File summary = "unmatched_filter_summary.json"
    }
}
