package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class PloidyTable {

    private final Map<String, PloidyRecord> samplePloidyMap;

    public PloidyTable(final Path path) {
        try {
            final TableReader<PloidyRecord> reader = TableUtils.reader(path, (columns, exceptionFactory) ->
                (dataLine) -> {
                    final String sample = dataLine.get(0);
                    final int numColumns = dataLine.columns().columnCount();
                    final List<String> contigs = columns.names().subList(1, numColumns);
                    final List<Integer> ploidy = IntStream.range(1, numColumns).map(i -> dataLine.getInt(i))
                            .boxed().collect(Collectors.toList());
                    return new PloidyRecord(sample, contigs, ploidy);
                }
            );
            samplePloidyMap = reader.stream().collect(Collectors.toMap(PloidyRecord::getSample, r -> r));
            reader.close();
        } catch (final IOException e) {
            throw new GATKException("IO error while reading ploidy table", e);
        }
    }

    public PloidyTable(final Map<String, Map<String, Integer>> samplePloidyMap) {
        this.samplePloidyMap = samplePloidyMap.entrySet().stream()
                .collect(Collectors.toMap(p -> p.getKey(), p -> new PloidyRecord(p.getKey(),
                        p.getValue().entrySet().stream().map(Map.Entry::getKey).collect(Collectors.toList()),
                        p.getValue().entrySet().stream().map(Map.Entry::getValue).collect(Collectors.toList()))));
    }

    public Integer get(final String sample, final String contig) {
        Utils.validateArg(samplePloidyMap.containsKey(sample), "Sample " + sample + " not found in ploidy records");
        return samplePloidyMap.get(sample).getPloidy(contig);
    }

    private static final class PloidyRecord {
        private final String sample;
        private final Map<String, Integer> ploidyMap;

        public PloidyRecord(final String sample, final List<String> contigs, final List<Integer> ploidy) {
            Utils.validateArg(contigs.size() == ploidy.size(), "Expected " + contigs.size() +
                    " ploidy values but found " + ploidy.size());
            this.sample = sample;
            this.ploidyMap = new HashMap<>();
            for (int i = 0; i < contigs.size(); i++) {
                ploidyMap.put(contigs.get(i), ploidy.get(i));
            }
        }

        public String getSample() {
            return sample;
        }

        public Integer getPloidy(final String contig) {
            Utils.validateArg(ploidyMap.containsKey(contig), "No ploidy entry for sample " + sample +
                    " at contig " + contig);
            return ploidyMap.get(contig);
        }
    }
}
