package org.broadinstitute.hellbender.tools.pon.allelic;

import com.google.common.collect.ImmutableMap;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountTableColumn;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Utility methods for reading/writing an {@link AllelicPanelOfNormals} from/to an {@link HDF5File}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class HDF5AllelicPoNUtils {
    private static final String ALLELIC_PANEL_GROUP_NAME = "/allelic_panel";
    private static final String GLOBAL_ALPHA_PATH = ALLELIC_PANEL_GROUP_NAME + "/global_alpha";
    private static final String GLOBAL_BETA_PATH = ALLELIC_PANEL_GROUP_NAME + "/global_beta";

    private static final String SNPS_PATH = ALLELIC_PANEL_GROUP_NAME + "/snps";

    private static final String NUM_SNP_COLUMNS_GROUP_NAME = "/num_snp_cols";
    private static final String NUM_SNP_COLUMNS_PATH = NUM_SNP_COLUMNS_GROUP_NAME + "/values";

    private static final String LOCAL_HYPERPARAMETERS_PATH = ALLELIC_PANEL_GROUP_NAME + "/local_alpha_beta";

    private static final int NUM_HYPERPARAMETERS = 2;

    private static final EnumMap<AllelicCountTableColumn, Integer> snpColumnToPoNIndex = new EnumMap<>(ImmutableMap.of(
            AllelicCountTableColumn.CONTIG, 0, AllelicCountTableColumn.POSITION, 1));
    private static final int NUM_SNP_COLUMNS = snpColumnToPoNIndex.size();

    private HDF5AllelicPoNUtils() {}

    /**
     * Reads an {@link AllelicPanelOfNormals} from an {@link HDF5File}.
     */
    static AllelicPanelOfNormals read(final HDF5File file) {
        if (file.isPresent(ALLELIC_PANEL_GROUP_NAME)) {
            //read global hyperparameters
            final double globalAlpha = file.readDouble(GLOBAL_ALPHA_PATH);
            final double globalBeta = file.readDouble(GLOBAL_BETA_PATH);
            final AllelicPanelOfNormals.HyperparameterValues globalHyperparameterValues =
                    new AllelicPanelOfNormals.HyperparameterValues(globalAlpha, globalBeta);
            //read SNPs
            final List<SimpleInterval> snps = readSNPs(file);
            //read local hyperparameters
            final List<AllelicPanelOfNormals.HyperparameterValues> localHyperparameterValues = readLocalhyperparameterValues(file);
            //check number of SNPs and local hyperparameters match
            if (snps.size() != localHyperparameterValues.size()) {
                throw new GATKException(String.format("Wrong number of elements in the SNPs and local hyperparameters recovered from file '%s': %d != %d", file.getFile(), snps.size(), localHyperparameterValues.size()));
            }
            //create site-to-hyperparameter map
            final Map<SimpleInterval, AllelicPanelOfNormals.HyperparameterValues> siteToHyperparameterPairMap = new HashMap<>();
            for (int i = 0; i < snps.size(); i++) {
                final SimpleInterval site = snps.get(i);
                final AllelicPanelOfNormals.HyperparameterValues hyperparameterValues = localHyperparameterValues.get(i);
                if (siteToHyperparameterPairMap.containsKey(site)) {
                    throw new UserException.BadInput("Allelic panel of normals contains duplicate sites.");
                } else {
                    siteToHyperparameterPairMap.put(site, hyperparameterValues);
                }
            }
            return new AllelicPanelOfNormals(globalHyperparameterValues, siteToHyperparameterPairMap);
        }
        //if HDF5 file does not contain path for allelic PoN, return empty PoN
        return AllelicPanelOfNormals.EMPTY_PON;
    }

    /**
     * Writes an {@link AllelicPanelOfNormals} to an {@link HDF5File}, which must be in appropriate {@link HDF5File.OpenMode}.
     */
    static void write(final HDF5File outputHDF5File, final AllelicPanelOfNormals allelicPoN) {
        Utils.nonNull(outputHDF5File);
        Utils.nonNull(allelicPoN);
        //write global hyperparameters
        final AllelicPanelOfNormals.HyperparameterValues globalHyperparameterValues = allelicPoN.getGlobalHyperparameterValues();
        final double globalAlpha = globalHyperparameterValues.getAlpha();
        final double globalBeta = globalHyperparameterValues.getBeta();
        outputHDF5File.makeDouble(GLOBAL_ALPHA_PATH, globalAlpha);
        outputHDF5File.makeDouble(GLOBAL_BETA_PATH, globalBeta);
        //write SNPs and local hyperparameters (both sorted by lexicographical contig order)
        final List<Map.Entry<SimpleInterval, AllelicPanelOfNormals.HyperparameterValues>> ponEntries = allelicPoN.getSortedMapEntries();
        final List<SimpleInterval> snps = ponEntries.stream().map(Map.Entry::getKey).collect(Collectors.toList());
        writeSNPValues(outputHDF5File, snps);
        final double[][] localHyperparameterMatrix = makeLocalHyperparameterMatrix(ponEntries);
        outputHDF5File.makeDoubleMatrix(LOCAL_HYPERPARAMETERS_PATH, localHyperparameterMatrix);
    }

    /**
     * Reads SNP contigs and positions from the appropriate path.
     */
    private static List<SimpleInterval> readSNPs(final HDF5File reader) {
        final String[][] values = reader.readStringMatrix(SNPS_PATH, NUM_SNP_COLUMNS_PATH);

        final int numSNPCols = (int) reader.readDouble(NUM_SNP_COLUMNS_PATH);
        final List<SimpleInterval> result = new ArrayList<>(values.length);
        for (final String[] row : values) {
            if (row.length != numSNPCols) {
                throw new GATKException(String.format("Wrong number of column elements in the SNPs recovered from file '%s': %d != %d", reader.getFile(), row.length, numSNPCols));
            }
            result.add(new SimpleInterval(row[0], Integer.parseInt(row[1]), Integer.parseInt(row[1])));
        }
        return result;
    }

    /**
     * Reads local hyperparameter values to the appropriate path.
     */
    private static List<AllelicPanelOfNormals.HyperparameterValues> readLocalhyperparameterValues(final HDF5File file) {
        final double[][] values = file.readDoubleMatrix(LOCAL_HYPERPARAMETERS_PATH);

        final List<AllelicPanelOfNormals.HyperparameterValues> result = new ArrayList<>(values.length);
        for (final double[] row : values) {
            if (row.length != NUM_HYPERPARAMETERS) {
                throw new GATKException(String.format("Wrong number of column elements in the local hyperparameters recovered from file '%s': %d != %d", file.getFile(), row.length, NUM_HYPERPARAMETERS));
            }
            result.add(new AllelicPanelOfNormals.HyperparameterValues(row[0], row[1]));
        }
        return result;
    }

    /**
     * Writes SNP contigs and positions to the appropriate path.
     */
    private static void writeSNPValues(final HDF5File file, final List<SimpleInterval> snps) {
        final String[][] snpValues = new String[snps.size()][NUM_SNP_COLUMNS];
        for (int i = 0; i < snps.size(); i++) {
            snpValues[i][snpColumnToPoNIndex.get(AllelicCountTableColumn.CONTIG)] = snps.get(i).getContig();
            snpValues[i][snpColumnToPoNIndex.get(AllelicCountTableColumn.POSITION)] = String.valueOf(snps.get(i).getStart());
        }
        file.makeStringMatrix(SNPS_PATH, snpValues, NUM_SNP_COLUMNS_PATH);
    }

    /**
     * Transforms map of local hyperparameters to a double matrix suitable for writing.
     */
    private static double[][] makeLocalHyperparameterMatrix(final List<Map.Entry<SimpleInterval, AllelicPanelOfNormals.HyperparameterValues>> ponEntries) {
        final double[][] localHyperparameterMatrix = new double[ponEntries.size()][NUM_HYPERPARAMETERS];
        for (int i = 0; i < ponEntries.size(); i++) {
            final AllelicPanelOfNormals.HyperparameterValues hyperparameterValues = ponEntries.get(i).getValue();
            localHyperparameterMatrix[i][0] = hyperparameterValues.getAlpha();
            localHyperparameterMatrix[i][1] = hyperparameterValues.getBeta();
        }
        return localHyperparameterMatrix;
    }
}
