package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.*;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

/**
 * Merge the stats output by scatters of a single Mutect2 job.
 */
@CommandLineProgramProperties(
        summary = "Merge the stats output by scatters of a single Mutect2 job",
        oneLineSummary = "Merge the stats output by scatters of a single Mutect2 job",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
public class MergeMutectStats extends CommandLineProgram {

    @Argument(fullName = Mutect2.MUTECT_STATS_SHORT_NAME, doc="Stats from Mutect2 scatters of a single tumor or tumor-normal pair")
    private Set<File> stats = new LinkedHashSet<>(0);

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="Output stats")
    private File outputStatsTable = null;

    public Object doWork() {
        final Map<String, ToDoubleFunction<List<Double>>> aggregations = new HashMap<>();
        aggregations.put(Mutect2Engine.CALLABLE_SITES_NAME, list -> list.stream().mapToDouble(x -> x).sum());

        final Map<String, List<Double>> scatteredStats = new HashMap<>();
        stats.stream().flatMap(file -> MutectStats.readFromFile(file).stream()).forEach(mutectStat -> {
            final String key = mutectStat.getStatistic();
            final double value = mutectStat.getValue();

            if (!scatteredStats.containsKey(key)) {
                scatteredStats.put(key, new ArrayList<>());
            }

            scatteredStats.get(key).add(value);
        });

        final List<MutectStats> aggregatedStats = scatteredStats.entrySet().stream().map(entry -> {
            final String stat = entry.getKey();
            final List<Double> values = entry.getValue();
            Utils.validate(aggregations.containsKey(stat), () -> "aggregations list missing key " + stat);
            final double aggregatedValue = aggregations.get(stat).applyAsDouble(values);
            return new MutectStats(stat, aggregatedValue);
        }).collect(Collectors.toList());

        MutectStats.writeToFile(aggregatedStats, outputStatsTable);

        return "SUCCESS";
    }
}

