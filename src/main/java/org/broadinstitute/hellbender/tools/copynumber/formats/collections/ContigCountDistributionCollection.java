package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.collections4.ListUtils;
import org.broadinstitute.hellbender.tools.copynumber.DetermineGermlineContigPloidy;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.Metadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.ContigCountDistribution;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.util.*;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Represents coverage distributions for each contig in an ordered set associated with a named sample.
 * Should only be used to write temporary files in {@link DetermineGermlineContigPloidy}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class ContigCountDistributionCollection extends AbstractRecordCollection<SampleLocatableMetadata, ContigCountDistribution> {
    private static final String CONTIG_TABLE_COLUMN = "CONTIG";

    private final int maximumCount;

    private static final Function<DataLine, ContigCountDistribution> CONTIG_COUNT_DISTRIBUTION_RECORD_FROM_DATA_LINE_DECODER = dataLine -> {
        final String contig = dataLine.get(CONTIG_TABLE_COLUMN);
        final Map<Integer, Integer> countDistribution = dataLine.columns().names().stream()
                .filter(c -> !c.equals(CONTIG_TABLE_COLUMN))
                .collect(Collectors.toMap(
                        Integer::parseInt,
                        dataLine::getInt));
        return new ContigCountDistribution(contig, countDistribution);
    };

    private static final BiConsumer<ContigCountDistribution, DataLine> CONTIG_COUNT_DISTRIBUTION_RECORD_TO_DATA_LINE_ENCODER = (contigCountDistribution, dataLine) -> {
        dataLine.append(contigCountDistribution.getContig());
        contigCountDistribution.getCountDistribution().forEach(
                (i, occurrences) -> dataLine.set(String.valueOf(i), occurrences)
        );
    };

    public ContigCountDistributionCollection(final SampleLocatableMetadata metadata,
                                             final List<ContigCountDistribution> contigCountDistributions) {
        super(metadata, contigCountDistributions, constructColumns(metadata, contigCountDistributions), CONTIG_COUNT_DISTRIBUTION_RECORD_FROM_DATA_LINE_DECODER, CONTIG_COUNT_DISTRIBUTION_RECORD_TO_DATA_LINE_ENCODER);
        maximumCount = contigCountDistributions.get(0).getCountDistribution().size() - 1;
    }

    private static TableColumnCollection constructColumns(final SampleLocatableMetadata metadata,
                                                          final List<ContigCountDistribution> contigCountDistributions) {
        Utils.nonNull(metadata);
        Utils.nonEmpty(contigCountDistributions);
        Utils.validateArg(contigCountDistributions.stream().map(ContigCountDistribution::getContig).distinct().count() == contigCountDistributions.size(),
                "List of contig count distributions cannot contain duplicate contigs.");
        final int maximumCount = contigCountDistributions.get(0).getCountDistribution().size() - 1;
        contigCountDistributions.forEach(
                ccd -> Utils.validateArg(ccd.getCountDistribution().size() == maximumCount + 1,
                        "Contig count distributions must all have the same maximum count."));
        return new TableColumnCollection(ListUtils.union(
                Collections.singletonList(CONTIG_TABLE_COLUMN),
                IntStream.range(0, maximumCount + 1).boxed().map(String::valueOf).collect(Collectors.toList())));
    }

    public ContigCountDistributionCollection(final SimpleCountCollection readCounts,
                                             final Set<SimpleInterval> intervalSubset,
                                             final int maximumCount) {
        this(Utils.nonNull(readCounts).getMetadata(), constructContigCountDistributions(readCounts, Utils.nonEmpty(intervalSubset), maximumCount));
    }

    private static List<ContigCountDistribution> constructContigCountDistributions(final SimpleCountCollection readCounts,
                                                                                   final Set<SimpleInterval> intervalSubset,
                                                                                   final int maximumCount) {
        Utils.nonNull(readCounts);
        Utils.nonEmpty(intervalSubset);
        ParamUtils.isPositiveOrZero(maximumCount, "Maximum count must be non-negative.");
        final Map<String, Map<Integer, Integer>> mapOfMaps = readCounts.getRecords().stream()
                .filter(c -> intervalSubset.contains(c.getInterval()))
                .collect(Collectors.groupingBy(
                        SimpleCount::getContig,
                        LinkedHashMap::new,
                        Collectors.groupingBy(
                                SimpleCount::getCount,
                                Collectors.summingInt(c -> 1))));
        return mapOfMaps.entrySet().stream()
                .map(e -> new ContigCountDistribution(
                        e.getKey(),
                        IntStream.range(0, maximumCount + 1).boxed()
                                .collect(Collectors.toMap(
                                        Function.identity(),
                                        i -> e.getValue().getOrDefault(i, 0)))))    //set the occurrences of unseen counts to zero
                .collect(Collectors.toList());
    }

    @Override
    protected Metadata.Type getMetadataType() {
        return Metadata.Type.SAMPLE_LOCATABLE;
    }

    public int getMaximumCount() {
        return maximumCount;
    }
}
