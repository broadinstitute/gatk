package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.tools.exome.HashedListTargetCollection;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.tools.walkers.mutect.FilterMutectCalls;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Given pileup data from {@link GetPileupSummaries}, calculates the fraction of reads coming from cross-sample contamination.
 *
 * <p>
 *     The resulting contamination table is used with {@link FilterMutectCalls}.
 * </p>
 *
 * <p>This tool and GetPileupSummaries together replace GATK3's ContEst.</p>
 *
 * <p>
 *     The resulting table provides the fraction contamination, one line per sample, e.g. SampleID--TAB--Contamination.
 *     The file has no header.
 * </p>
 *
 * <h3>Example</h3>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" CalculateContamination \
 *   -I pileups.table \
 *   -O contamination.table
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Calculate contamination",
        oneLineSummary = "Calculate contamination",
        programGroup = VariantProgramGroup.class
)
@DocumentedFeature
public class CalculateContamination extends CommandLineProgram {

    private static final Logger logger = LogManager.getLogger(CalculateContamination.class);

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc="The input table", optional = false)
    private File inputPileupSummariesTable;

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output table", optional = false)
    private final File outputTable = null;

    private static final int CNV_SCALE = 1000000;

    @Override
    public Object doWork() {
        final List<PileupSummary> pileupSummaries = PileupSummary.readPileupSummaries(inputPileupSummariesTable);
        final List<PileupSummary> homAltSites = findConfidentHomAltSites(pileupSummaries);
        final Pair<Double, Double> contaminationAndError = homAltSites.isEmpty() ? Pair.of(0.0, 0.0) : calculateContamination(homAltSites);
        ContaminationRecord.writeContaminationTable(Arrays.asList(new ContaminationRecord(ContaminationRecord.Level.WHOLE_BAM.toString(), contaminationAndError.getLeft(), contaminationAndError.getRight())), outputTable);

        return "SUCCESS";
    }

    private Pair<Double, Double> calculateContamination(List<PileupSummary> homAltSites) {

        final long totalReadCount = homAltSites.stream().mapToLong(PileupSummary::getTotalCount).sum();
        final long totalRefCount = homAltSites.stream().mapToLong(PileupSummary::getRefCount).sum();

        // if eg ref is A, alt is C, then # of ref reads due to error is roughly (# of G read + # of T reads)/2
        final long errorRefCount = homAltSites.stream().mapToLong(PileupSummary::getOtherAltCount).sum() / 2;
        final long contaminationRefCount = Math.max(totalRefCount - errorRefCount, 0);
        final double totalDepthWeightedByRefFraction = homAltSites.stream()
                .mapToDouble(ps -> ps.getTotalCount() * (1 - ps.getAlleleFrequency()))
                .sum();
        final double contamination = contaminationRefCount / totalDepthWeightedByRefFraction;
        final double standardError = Math.sqrt(contamination / totalDepthWeightedByRefFraction);

        logger.info(String.format("In %d homozygous variant sites we find %d reference reads due to contamination and %d" +
                        " due to to sequencing error out of a total %d reads.", homAltSites.size(), contaminationRefCount, errorRefCount, totalReadCount));
        logger.info(String.format("Based on population data, we would expect %d reference reads in a contaminant with the same depths at these sites.", (long) totalDepthWeightedByRefFraction));
        logger.info(String.format("Therefore, we estimate a contamination of %.3f.", contamination));
        logger.info(String.format("The error bars on this estimate are %.5f.", standardError));
        return Pair.of(contamination, standardError);
    }

    private static List<PileupSummary> findConfidentHomAltSites(List<PileupSummary> sites) {
        if (sites.isEmpty()) {
            return new ArrayList<>();
        }
        final TargetCollection<PileupSummary> tc = new HashedListTargetCollection<>(sites);
        final double averageCoverage = sites.stream().mapToInt(PileupSummary::getTotalCount).average().getAsDouble();
        final List<Double> smoothedCopyRatios = new ArrayList<>();
        final List<Double> hetRatios = new ArrayList<>();

        for (final PileupSummary site : sites) {
            final SimpleInterval nearbySpan = new SimpleInterval(site.getContig(), Math.max(1, site.getStart() - CNV_SCALE), site.getEnd() + CNV_SCALE);
            final List<PileupSummary> nearbySites = tc.targets(nearbySpan);

            final double averageNearbyCopyRatio = nearbySites.stream().mapToDouble(s -> s.getTotalCount()/averageCoverage).average().orElseGet(() -> 0);
            smoothedCopyRatios.add(averageNearbyCopyRatio);

            final double expectedNumberOfNearbyHets = nearbySites.stream().mapToDouble(PileupSummary::getAlleleFrequency).map(x -> 2*x*(1-x)).sum();
            final long numberOfNearbyHets = nearbySites.stream().mapToDouble(PileupSummary::getAltFraction).filter(x -> 0.4 < x && x < 0.6).count();
            final double hetRatio = numberOfNearbyHets / expectedNumberOfNearbyHets;
            hetRatios.add(hetRatio);
        }

        final double medianSmoothedCopyRatio = new Median().evaluate(smoothedCopyRatios.stream().mapToDouble(x->x).toArray());
        final List<Integer> indicesWithAnomalousCopyRatio = IntStream.range(0, sites.size())
                .filter(n -> smoothedCopyRatios.get(n) < 0.8 * medianSmoothedCopyRatio || smoothedCopyRatios.get(n) > 2 *medianSmoothedCopyRatio)
                .boxed().collect(Collectors.toList());

        final double meanHetRatio = hetRatios.stream().mapToDouble(x->x).average().getAsDouble();
        final List<Integer> indicesWithLossOfHeterozygosity = IntStream.range(0, sites.size())
                .filter(n -> hetRatios.get(n) < meanHetRatio * 0.5)
                .boxed().collect(Collectors.toList());

        //TODO: as extra security, filter out sites that are near too many hom alts

        logger.info(String.format("Excluding %d sites with low or high copy ratio and %d sites with potential loss of heterozygosity",
                indicesWithAnomalousCopyRatio.size(), indicesWithLossOfHeterozygosity.size()));

        logger.info(String.format("The average ratio of hets within distance %d to theoretically expected number of hets is %.3f", CNV_SCALE, meanHetRatio));

        final Set<Integer> badSites = new TreeSet<>();
        badSites.addAll(indicesWithAnomalousCopyRatio);
        badSites.addAll(indicesWithLossOfHeterozygosity);

        return IntStream.range(0, sites.size())
                .filter(n -> !badSites.contains(n))
                .mapToObj(sites::get)
                .filter(s -> s.getAltFraction() > 0.8)
                .collect(Collectors.toList());
    }
}
