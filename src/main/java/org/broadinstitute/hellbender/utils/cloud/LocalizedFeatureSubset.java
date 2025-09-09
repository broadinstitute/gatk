package org.broadinstitute.hellbender.utils.cloud;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.sv.SVFeature;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.FeatureOutputCodec;
import org.broadinstitute.hellbender.utils.codecs.FeatureOutputCodecFinder;
import org.broadinstitute.hellbender.utils.codecs.FeatureSink;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Spliterators;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

public class LocalizedFeatureSubset<T extends SVFeature> {

    private FeatureDataSource<T> localizedSource;
    private final long recordCount;

    @SuppressWarnings("unchecked")
    public LocalizedFeatureSubset(final GATKPath featurePath, final Collection<SimpleInterval> intervals,
                                  final String name, final int localQueryLookaheadBases, final int cloudPrefetchBuffer,
                                  final int cloudIndexPrefetchBuffer, final List<String> samples,
                                  final int compressionLevel, final SAMSequenceDictionary dictionary) {
        final FeatureOutputCodec<? extends Feature, ? extends FeatureSink<? extends Feature>> outputCodec = FeatureOutputCodecFinder.find(featurePath);
        final Class<? extends Feature> outputClass = outputCodec.getFeatureType();
        if ( !SVFeature.class.isAssignableFrom(outputClass) ) {
            throw new UserException("F " + featurePath + " implies Feature subtype " +
                    outputClass.getSimpleName() + " but this class requires an SVFeature subtype.");
        }
        // Just use input filename as the extension to avoid parsing issues
        final String extension = featurePath.toPath().getFileName().toString();
        final File localizedFeaturesFile = IOUtils.createTempFile(name + "_cache", extension);

        // the validity of this cast was checked in the constructor
        try (final FeatureDataSource<T> originalSource = new FeatureDataSource<T>(featurePath.toString(), name, 0, outputClass, cloudPrefetchBuffer, cloudIndexPrefetchBuffer);
            final FeatureSink<SVFeature> outputSink = (FeatureSink<SVFeature>) outputCodec.makeSink(new GATKPath(localizedFeaturesFile.getPath()), dictionary, samples, compressionLevel)) {
            final Iterator<T> iter = IntervalUtils.sortAndMergeIntervalsToStream(intervals, dictionary, IntervalMergingRule.ALL)
                    .map(originalSource::query)
                    .flatMap(this::iteratorToStream)
                    .iterator();
            long count = 0;
            while (iter.hasNext()) {
                outputSink.write(iter.next());
                count++;
            }
            recordCount = count;
        }
        localizedSource = new FeatureDataSource<>(localizedFeaturesFile, name + "_local", localQueryLookaheadBases);
    }

    public long getRecordCount() {
        return recordCount;
    }

    public FeatureDataSource<T> getLocalizedDataSource() {
        return localizedSource;
    }

    public Iterator<T> query(final SimpleInterval interval) {
        return localizedSource.query(interval);
    }

    private Stream<T> iteratorToStream(final Iterator<T> iter) {
        return StreamSupport.stream(Spliterators.spliteratorUnknownSize(iter, 0), false);
    }

}
