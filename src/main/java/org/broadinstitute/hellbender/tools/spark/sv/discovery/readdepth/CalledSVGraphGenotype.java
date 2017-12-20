package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Called haplotypes class with associated list of coordinate-defined paths (for output)
 */
public final class CalledSVGraphGenotype extends SVGraphGenotype {

    private final List<CoordinateSVGraphPath> coordinatePaths;

    public final static String GRAPH_ID_COLUMN = "GID";
    public final static String GENOTYPE_ID_COLUMN = "GTID";
    public final static String HAPLOTYPE_ID_COLUMN = "HID";
    public final static String LIKELIHOOD_COLUMN = "LIK";
    public final static String PROBABILITY_COLUMN = "P";
    public final static String PATHS_COLUMN = "PATHS";

    public CalledSVGraphGenotype(final SVGraphGenotype genotype, final SVGraph graph) {
        super(genotype.getGroupId(), genotype.getGenotypeId(), genotype.getLikelihood());
        this.probability = genotype.getProbability();
        this.coordinatePaths = genotype.getHaplotypes().stream().map(path -> path.convertToCoordinatePath(graph)).collect(Collectors.toList());
    }

    public static String bedHeader() {
        return String.join("\t", Arrays.asList(GRAPH_ID_COLUMN, GENOTYPE_ID_COLUMN, HAPLOTYPE_ID_COLUMN, LIKELIHOOD_COLUMN, PROBABILITY_COLUMN, PATHS_COLUMN));
    }

    public List<String> bedStrings() {
        final List<String> pathStrings = coordinatePaths.stream().map(CoordinateSVGraphPath::toString).collect(Collectors.toList());
        final List<String> lines = new ArrayList<>(pathStrings.size());
        for (int haplotypeId = 0; haplotypeId < pathStrings.size(); haplotypeId++) {
            lines.add(getGroupId() + "\t" + getGenotypeId() + "\t" + haplotypeId + "\t" +
                    getLikelihood() + "\t" + getProbability() + "\t" + pathStrings.get(haplotypeId));
        }
        return lines;
    }

}
