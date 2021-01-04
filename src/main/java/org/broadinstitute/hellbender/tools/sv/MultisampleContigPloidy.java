package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledContigPloidyCollection;
import org.broadinstitute.hellbender.tools.walkers.sv.SVModelToolUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

public class MultisampleContigPloidy {

    private final Map<SampleContigPair, Integer> contigPloidyMap;

    public MultisampleContigPloidy(final List<GATKPath> contigPloidyCallFilePaths) {
        Utils.nonNull(contigPloidyCallFilePaths);
        final Collection<CalledContigPloidyCollection> contigPloidyCollections = contigPloidyCallFilePaths.stream()
                .map(p -> new CalledContigPloidyCollection(p.toPath().toFile()))
                .collect(Collectors.toList());
        final List<String> samples = contigPloidyCollections.stream().map(CalledContigPloidyCollection::getMetadata)
                .map(m -> m.getSampleName()).collect(Collectors.toList());
        SVModelToolUtils.assertDistinctSampleIds(samples);
        final Map<String, CalledContigPloidyCollection> sampleToCollectionMap = contigPloidyCollections.stream()
                .collect(Collectors.toMap(c -> c.getMetadata().getSampleName(), p -> p));
        contigPloidyMap = getInternalMap(sampleToCollectionMap);
    }

    public MultisampleContigPloidy(final Map<String, CalledContigPloidyCollection> sampleToCollectionMap) {
        Utils.nonNull(sampleToCollectionMap);
        contigPloidyMap = getInternalMap(sampleToCollectionMap);
    }

    private Map<SampleContigPair, Integer> getInternalMap(final Map<String, CalledContigPloidyCollection> contigPloidyCollectionMap) {
        return contigPloidyCollectionMap.entrySet().stream()
                .flatMap(e -> e.getValue().getRecords().stream().map(r -> new AbstractMap.SimpleImmutableEntry<>(new SampleContigPair(e.getKey(), r.getContig()), r.getPloidy())))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
    }

    public Set<String> getSampleSet() {
        return contigPloidyMap.keySet().stream().map(SampleContigPair::getSample).collect(Collectors.toSet());
    }

    public Set<String> getContigSet() {
        return contigPloidyMap.keySet().stream().map(SampleContigPair::getSample).collect(Collectors.toSet());
    }

    public Integer get(final String sample, final String contig) {
        return contigPloidyMap.get(new SampleContigPair(sample, contig));
    }

    private final static class SampleContigPair {
        private final String sample;
        private final String contig;
        public SampleContigPair(final String sample, final String contig) {
            this.sample = sample;
            this.contig = contig;
        }

        public String getSample() {
            return sample;
        }

        public String getContig() {
            return contig;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof SampleContigPair)) return false;
            SampleContigPair that = (SampleContigPair) o;
            return Objects.equals(sample, that.sample) && Objects.equals(contig, that.contig);
        }

        @Override
        public int hashCode() {
            return Objects.hash(sample, contig);
        }
    }
}
