package org.broadinstitute.hellbender.tools.spark.sv;

import avro.shaded.com.google.common.collect.Sets;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedAssembly;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;

/**
 * TODO this is just a brute sh_playground_ctx branch copy&paste... need to change it to make
 * reference to steves code once gets merged into master.
 */
public class PlaygroundExtract {


    public static List<Tuple2<List<AlignmentInterval>, Double>> scoreAllConfigurations(final AlignedContig alignedContig,
                                                                                       final Set<String> canonicalChromosomes,
                                                                                       final int maxCanonicalChrAlignerScore) {

        // speed up if number of alignments is too high (>10)
        // if mapped to canonical chromosomes, MQ must be >10; otherwise, must have AS higher than max canonical aligner score
        final List<AlignmentInterval> alignmentIntervals;
        if (alignedContig.alignmentIntervals.size() > 10) {
            alignmentIntervals = alignedContig.alignmentIntervals.stream()
                    .filter(alignmentInterval -> (!canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig())
                            && alignmentInterval.alnScore > maxCanonicalChrAlignerScore)
                            || alignmentInterval.mapQual>10)
                    .collect(Collectors.toList());
        } else {
            alignmentIntervals = alignedContig.alignmentIntervals;
        }

        final int newMaxCanonicalChrAlignerScore = alignmentIntervals.stream()
                .filter(alignmentInterval -> canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig()))
                .mapToInt(ai -> ai.alnScore).max().orElse(0); // possible that no mapping to canonical chromosomes

        final Comparator<AlignmentInterval> alignedIntervalComparator =
                Comparator.comparingInt(x -> x.referenceSpan.getStart());
        // group 2: generate, and score configurations
        final List<List<AlignmentInterval>> allConfigurations = Sets.powerSet(new HashSet<>(alignmentIntervals))
                .stream().map(ArrayList::new)
                // make sure within each configuration, alignments would be sorted as they would be in a corresponding AlignedContig
                .map(ls -> ls.stream().sorted(alignedIntervalComparator).collect(Collectors.toList()))
                .collect(Collectors.toList());

        return allConfigurations.stream()
                .map(configuration ->
                        new Tuple2<>(configuration, computeScoreOfConfiguration(configuration, canonicalChromosomes, newMaxCanonicalChrAlignerScore)))
                .collect(Collectors.toCollection(() -> new ArrayList<>(allConfigurations.size())));
    }

    public static Tuple2<List<AlignmentInterval>, Double> pickAnyBestConfiguration(final AlignedContig alignedContig, final Set<String> canonicalChromosomes) {
        final int maxCanonicalChrAlignerScore = alignedContig.alignmentIntervals.stream()
                .filter(alignmentInterval -> canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig()))
                .mapToInt(ai -> ai.alnScore).max().orElse(0); // possible that no mapping to canonical chromosomes


        final List<Tuple2<List<AlignmentInterval>, Double>> scoredConfiguration = scoreAllConfigurations(alignedContig, canonicalChromosomes, maxCanonicalChrAlignerScore);
        return scoredConfiguration.stream()
                .sorted(Comparator.comparingDouble(t -> -t._2().doubleValue()))
                .findFirst().orElse(new Tuple2<>(Collections.emptyList(), Double.NEGATIVE_INFINITY));
    }

    static List<List<AlignmentInterval>> pickBestConfigurations(final AlignedContig alignedContig,
                                                                                final Set<String> canonicalChromosomes) {

        // group 1: get max aligner score of mappings to canonical chromosomes and speed up in case of too many mappings
        final int maxCanonicalChrAlignerScore = alignedContig.alignmentIntervals.stream()
                .filter(alignmentInterval -> canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig()))
                .mapToInt(ai -> ai.alnScore).max().orElse(0); // possible that no mapping to canonical chromosomes

        final  List<Tuple2<List<AlignmentInterval>, Double>> scoredConfiguration = scoreAllConfigurations(alignedContig, canonicalChromosomes, maxCanonicalChrAlignerScore);

        // group 3: pick the best-scored configuration(s) (if multiple configuration equally-good scored, return all of them)
        final double maxScore = scoredConfiguration.stream().mapToDouble(Tuple2::_2).max()
                .orElseThrow(() -> new GATKException("Cannot find best-scoring configuration on alignments of contig: " + alignedContig.contigName));

        return scoredConfiguration.stream()
                .filter(t -> t._2() >= maxScore)
                .map(Tuple2::_1)
                .collect(Collectors.toList());
    }

    public static double computeScoreOfConfiguration(final List<AlignmentInterval> configuration,
                                              final Set<String> canonicalChromosomes,
                                              final int maxCanonicalChrAlignerScore) {

        double tigExplainQual =
                configuration.stream()
                        .mapToDouble(alignmentInterval -> {
                            final int len = alignmentInterval.endInAssembledContig - alignmentInterval.startInAssembledContig + 1;
                            final double weight;
                            if (canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig())) {
                                weight = alignmentInterval.mapQual/60.0;
                            } else {
                                weight = Math.max(alignmentInterval.mapQual/60.0,
                                        alignmentInterval.alnScore > maxCanonicalChrAlignerScore ? 1 : 0);
                            }
                            return weight * len;
                        })
                        .sum();

        int redundancy = 0;
        for (int i = 0; i < configuration.size() -1 ; ++i) {
            final int b1 = configuration.get(i).endInAssembledContig;
            for (int j = i + 1; j < configuration.size(); ++j) {
                final int a2 = configuration.get(j).startInAssembledContig;
                final int overlap = Math.max(b1 - a2 + 1, 0);
                redundancy += overlap;
            }
        }

        return tigExplainQual - redundancy;
    }

}
