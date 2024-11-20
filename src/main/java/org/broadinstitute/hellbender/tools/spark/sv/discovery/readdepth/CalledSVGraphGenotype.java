package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Called haplotypes class with associated list of coordinate-defined paths (for output)
 */
public final class CalledSVGraphGenotype extends SVGraphGenotype {

    private final List<IndexedSVGraphPath> coordinatePaths;
    private final SVGraph graph;

    public final static String GRAPH_ID_COLUMN = "GID";
    public final static String GENOTYPE_ID_COLUMN = "GTID";
    public final static String HAPLOTYPE_ID_COLUMN = "HID";
    public final static String DEPTH_LIKELIHOOD_COLUMN = "DL";
    public final static String DEPTH_PROBABILITY_COLUMN = "DP";
    public final static String EVIDENCE_PROBABILITY_COLUMN = "EP";
    public final static String PROBABILITY_COLUMN = "P";
    public final static String PATHS_COLUMN = "PATHS";

    public CalledSVGraphGenotype(final SVGraphGenotype genotype, final SVGraph graph) {
        super(genotype.getGroupId(), genotype.getGenotypeId(), genotype.getHaplotypes());
        this.probability = genotype.getProbability();
        this.depthLikelihood = genotype.getDepthLikelihood();
        this.depthProbability = genotype.getDepthProbability();
        this.evidenceProbability = genotype.getEvidenceProbability();
        this.coordinatePaths = genotype.getHaplotypes();
        this.graph = graph;
    }

    public static String bedHeader() {
        return String.join("\t", Arrays.asList(GRAPH_ID_COLUMN, GENOTYPE_ID_COLUMN, HAPLOTYPE_ID_COLUMN, PROBABILITY_COLUMN, DEPTH_LIKELIHOOD_COLUMN, DEPTH_PROBABILITY_COLUMN, EVIDENCE_PROBABILITY_COLUMN, PATHS_COLUMN));
    }

    public List<String> bedStrings() {
        final List<String> pathStrings = coordinatePaths.stream().map(path -> path.convertToCoordinatePath(graph).toString()).collect(Collectors.toList());
        final List<String> lines = new ArrayList<>(pathStrings.size());
        for (int haplotypeId = 0; haplotypeId < pathStrings.size(); haplotypeId++) {
            lines.add(getGroupId() + "\t" + getGenotypeId() + "\t" + haplotypeId + "\t" + getProbability() + "\t" +
                    getDepthLikelihood() + "\t" + getDepthProbability() + "\t" + getEvidenceProbability() + "\t" + pathStrings.get(haplotypeId));
        }
        return lines;
    }

}
