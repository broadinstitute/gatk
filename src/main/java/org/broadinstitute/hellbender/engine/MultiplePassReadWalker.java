package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Use this extension of a ReadWalker when you need to process the reads multiple times.
 * Implement the method traverseReads(), calling forEachRead(func) as many times as you need to.
 * The functional object you provide to forEachRead will be presented with each read.
 */
public abstract class MultiplePassReadWalker extends ReadWalker {
    private SAMSequenceDictionary dictionary;
    private CountingReadFilter countedFilter;
    private SamReaderFactory samReaderFactory;

    @FunctionalInterface
    public interface GATKApply {
        void apply( GATKRead read, ReferenceContext reference, FeatureContext features );
    }

    public abstract void traverseReads();

    public void forEachRead( final GATKApply func ) {
        try ( final ReadsDataSource readsSource =
                      new ReadsDataSource(
                              readArguments.getReadPaths(),
                              readArguments.getReadIndexPaths(),
                              samReaderFactory,
                              cloudPrefetchBuffer,
                              (cloudIndexPrefetchBuffer < 0 ? cloudPrefetchBuffer : cloudIndexPrefetchBuffer)) ) {

            if ( hasUserSuppliedIntervals() ) {
                readsSource.setTraversalBounds(intervalArgumentCollection.getTraversalParameters(dictionary));
            }

            final ReadTransformer preTransformer = makePreReadFilterTransformer();
            final ReadTransformer postTransformer = makePostReadFilterTransformer();
            Utils.stream(readsSource).map(preTransformer).filter(countedFilter).map(postTransformer).forEach( read -> {
                final SimpleInterval readInterval = getReadInterval(read);
                func.apply(
                        read,
                        new ReferenceContext(reference, readInterval), // will be empty if reference or readInterval is null
                        new FeatureContext(features, readInterval));   // will be empty if features or readInterval is null
                progressMeter.update(readInterval); });
        }
    }

    @Override
    public final void traverse() {
        dictionary = reads.getHeader().getSequenceDictionary();
        reads.close();
        reads = null;

        countedFilter = makeReadFilter();
        SamReaderFactory factory =
                SamReaderFactory.makeDefault().validationStringency(readArguments.getReadValidationStringency());
        if ( hasReference() ) { // pass in reference if available, because CRAM files need it
            factory = factory.referenceSequence(referenceArguments.getReferencePath());
        }
        if ( intervalArgumentCollection.intervalsSpecified() && !disableBamIndexCaching ) {
            factory = factory.enable(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES);
        }
        samReaderFactory = factory;
        traverseReads();
        logger.info(countedFilter.getSummaryLine());
    }

    // not used.  here because we extend ReadWalker, which requires it.
    @Override
    public final void apply( final GATKRead read,
                             final ReferenceContext referenceContext,
                             final FeatureContext featureContext ) {
        throw new GATKException("apply called unexpectedly");
    }
}
