package org.broadinstitute.hellbender.tools.walkers.sv;

import avro.shaded.com.google.common.collect.Sets;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public final class SVModelToolUtils {

    public static Map<String, Double> readSampleDepthTable(final GATKPath path) {
        try {
            final List<Map.Entry<String, Double>> entries = new BufferedReader(new FileReader(path.toPath().toFile()))
                    .lines()
                    .map(s -> s.split("\t"))
                    .map(arr -> new AbstractMap.SimpleImmutableEntry<>(arr[0], Double.valueOf(arr[1])))
                    .collect(Collectors.toList());
            final List<String> samples = entries.stream().map(Map.Entry::getKey).collect(Collectors.toList());
            assertDistinctSampleIds(samples);
            return entries.stream().collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
        } catch (final IOException e) {
            throw new GATKException("Error reading sample depth file", e);
        }
    }

    public static Set<String> assertDistinctSampleIds(final Collection<String> samples) {
        final Set<String> set = new HashSet<>(samples);
        if (set.size() != samples.size()) {
            throw new UserException.BadInput("Sample ids are not distinct");
        }
        return set;
    }

    public static List<String> negotiateSampleSets(final List<Set<String>> sampleSets, final Logger logger) {
        Utils.nonNull(sampleSets);
        Utils.nonEmpty(sampleSets);
        Utils.containsNoNull(sampleSets, "Null sample set");
        sampleSets.stream().forEach(Utils::nonEmpty);
        final List<String> setSizes = sampleSets.stream().map(Set::size).map(String::valueOf).collect(Collectors.toList());
        final String inputSizesString = String.join(", ", setSizes);
        Set<String> resultSet = new HashSet<>(sampleSets.iterator().next());
        for (final Set<String> set : sampleSets) {
            resultSet = Sets.intersection(resultSet, set);
        }
        if (resultSet.isEmpty()) {
            throw new UserException.BadInput("Inputs had no samples in common");
        }
        logger.info("Using sample set of size " + resultSet.size() + " from input(s) of size " + inputSizesString);
        return resultSet.stream().sorted().collect(Collectors.toList());
    }

    public static File createSampleDepthTable(final List<String> orderedSampleList, final Map<String, Double> sampleDepthMap) {
        final List<String> depthTableFileLines = orderedSampleList.stream().map(s -> s + "\t" + sampleDepthMap.get(s)).collect(Collectors.toList());
        return IOUtils.writeTempFile(depthTableFileLines, "model.sample_depth_table", ".tmp");
    }

    public static File createSampleFile(final List<String> orderedSampleList) {
        return IOUtils.writeTempFile(orderedSampleList, "model.samples", ".tmp");
    }

}
