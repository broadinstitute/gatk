package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.avro.generic.GenericRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.variantdb.SchemaUtils;
import org.broadinstitute.hellbender.utils.bigquery.StorageAPIAvroReader;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

public class FilterSensitivityTools {
    private static final Logger logger = LogManager.getLogger(FilterSensitivityTools.class);

    private static double DEFAULT_TRUTH_SENSITIVITY_THRESHOLD_SNPS = 99.7;
    private static double DEFAULT_TRUTH_SENSITIVITY_THRESHOLD_INDELS = 99.0;

    public static void validateFilteringCutoffs(
            Double truthSensitivitySNPThreshold,
            Double truthSensitivityINDELThreshold,
            Double vqsLodSNPThreshold,
            Double vqsLodINDELThreshold,
            String tranchesTableName) {
        boolean snpTruthSensIsDefined = truthSensitivitySNPThreshold != null;
        boolean indelTruthSensIsDefined = truthSensitivityINDELThreshold != null;

        boolean snpVqslodIsDefined = vqsLodSNPThreshold != null;
        boolean indelVqslodIsDefined = vqsLodINDELThreshold != null;

        if (snpTruthSensIsDefined ^ indelTruthSensIsDefined) {
            throw new UserException("If one of (--snps-truth-sensitivity-filter-level, --indels-truth-sensitivity-filter-level) is provided, both must be provided.");
        } else if (snpTruthSensIsDefined) {
            // if the user specifies both truth sensitivity thresholds and lod cutoffs then throw a user error
            if (snpVqslodIsDefined || indelVqslodIsDefined) {
                throw new UserException("Arguments --[snps/indels]-truth-sensitivity-filter-level and --[snps/indels]-lod-score-cutoff are mutually exclusive. Please only specify one set of options.");
            }
        } else if (snpVqslodIsDefined ^ indelVqslodIsDefined) {
            throw new UserException("If one of (--snps-lod-score-cutoff, --indels-lod-score-cutoff) is provided, both must be provided.");
        }

        if (!snpVqslodIsDefined) {
            // we will need to use tranches to look up vqslod thresholds; therefore user must supply
            // a tranches table
            if (tranchesTableName == null) {
                throw new UserException("Unless using lod score cutoffs (advanced), you must provide a tranches table using the argument --tranches-table.");
            }
        }
    }

    public static Map<String, Map<Double, Double>> getTrancheMaps(String filterSetName, String tranchesTableName, String projectID) {
        // get tranches from BigQuery
        final String restrictionWithFilterSetName = SchemaUtils.FILTER_SET_NAME + " = '" + filterSetName + "'";
        final TableReference tranchesTableRef = new TableReference(tranchesTableName, SchemaUtils.TRANCHE_FIELDS);
        final StorageAPIAvroReader tranchesTableAvroReader = new StorageAPIAvroReader(tranchesTableRef, restrictionWithFilterSetName, projectID);

        Map<Double, Double> snpTranches = new TreeMap<>();
        Map<Double, Double> indelTranches = new TreeMap<>();

        for ( final GenericRecord queryRow : tranchesTableAvroReader ) {
            double targetTruthSensitivity = Double.parseDouble(queryRow.get(SchemaUtils.TARGET_TRUTH_SENSITIVITY).toString());
            double minVqslod = Double.parseDouble(queryRow.get(SchemaUtils.MIN_VQSLOD).toString());
            String model = queryRow.get(SchemaUtils.TRANCHE_MODEL).toString();
            if (GATKVCFConstants.SNP.equals(model)) {
                snpTranches.put(targetTruthSensitivity, minVqslod);
            } else if (GATKVCFConstants.INDEL.equals(model)) {
                indelTranches.put(targetTruthSensitivity, minVqslod);
            } else {
                // should this be an exception instead?
                logger.warn("Ignoring unrecognized model '" + model + "' found in tranches table " + tranchesTableName + ", filter set " + filterSetName);
            }
        }

        tranchesTableAvroReader.close();

        Map<String, Map<Double, Double>> trancheMaps = new HashMap<>();
        trancheMaps.put(GATKVCFConstants.SNP, snpTranches);
        trancheMaps.put(GATKVCFConstants.INDEL, indelTranches);

        return trancheMaps;
    }


    public static Double getVqslodThreshold(Map<Double, Double> trancheMap, Double truthSensitivityThreshold, String model) {

        if (truthSensitivityThreshold == null) {  // at this point, we know that all vqsr threshold inputs are null, so use defaults
            truthSensitivityThreshold = GATKVCFConstants.SNP.contains(model) ? DEFAULT_TRUTH_SENSITIVITY_THRESHOLD_SNPS : DEFAULT_TRUTH_SENSITIVITY_THRESHOLD_INDELS;
            logger.info("No filter thresholds supplied; using default " + model + " truth sensitivity threshold of " + truthSensitivityThreshold);
        }

        logger.info("Retrieving the min vqslod threshold from tranches for " + model + "s with truth sensitivity threshold of " + truthSensitivityThreshold);
        Double vqslodThreshold = getVqslodThresholdFromTranches(trancheMap, truthSensitivityThreshold);
        logger.info("Passing all " + model + " variants with VQSLOD >= " + vqslodThreshold);
        return vqslodThreshold;
    }


    public static Double getVqslodThresholdFromTranches(Map<Double, Double> trancheMap, Double truthSensitivityThreshold) {
        // We want to find the tranche with the smallest target_truth_sensitivity that is
        // equal to or greater than our truthSensitivityThreshold.
        // e.g. if truthSensitivitySNPThreshold is 99.8 and we have tranches with target_truth_sensitivities
        // of 99.5, 99.7, 99.9, and 100.0, we want the 99.9 sensitivity tranche.

        Double effectiveSensitivity = trancheMap.keySet().stream().filter(ts -> ts >= truthSensitivityThreshold).findFirst().orElse(null);

        if (effectiveSensitivity == null) {
            throw new UserException("No tranches found with target_truth_sensitivity >= " + truthSensitivityThreshold);
        }

        Double minVqslod = trancheMap.get(effectiveSensitivity);

        logger.info("Using tranche defined by sensitivity " + effectiveSensitivity + " and VQSLOD >= " + minVqslod + "; keeping all variants in this tranche.");

        return minVqslod;
    }


    public static VCFHeaderLine getVqsLodHeader(Double vqsLodThreshold, String model) {
        return new VCFFilterHeaderLine(GATKVCFConstants.VQSR_FAILURE_PREFIX + model,
                "Site failed " + model + " model VQSLOD cutoff of " + vqsLodThreshold.toString());
    }

    public static VCFHeaderLine getTruthSensitivityHeader(Double truthSensitivityThreshold, Double vqsLodThreshold, String model) {
        if (truthSensitivityThreshold == null) {  // at this point, we know that all vqsr threshold inputs are null, so use defaults
            truthSensitivityThreshold = GATKVCFConstants.SNP.contains(model) ? DEFAULT_TRUTH_SENSITIVITY_THRESHOLD_SNPS : DEFAULT_TRUTH_SENSITIVITY_THRESHOLD_INDELS;
        }
        return new VCFFilterHeaderLine(GATKVCFConstants.VQSR_FAILURE_PREFIX + model,
                "Site failed " + model + " model sensitivity cutoff (" + truthSensitivityThreshold.toString() + "), corresponding with VQSLOD cutoff of " + vqsLodThreshold.toString());
    }
}
