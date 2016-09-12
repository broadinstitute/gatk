package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.samtools.util.LocationAware;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.index.tabix.TabixFormat;
import org.broadinstitute.hellbender.engine.ProgressMeter;

import java.io.IOException;
import java.io.InputStream;

/**
 * This class is useful when we want to report progress when indexing. The indexing code is in htsjdk and not available to us directly.
 * The workaround is to make a special decorator 'codec' that gets called on every decoded Feature and this can use used to track progress.
 * The codec delegates all calls and reports progress as it goes.
 */
public final class ProgressReportingDelegatingCodec<A extends Feature, B> implements FeatureCodec<A, B> {
    private final FeatureCodec<A, B> delegatee;

    private final ProgressMeter pm;

    //Note: this default constructor is needed for the FeatureManager when it loads codecs.
    @SuppressWarnings("unused")
    public ProgressReportingDelegatingCodec(){
        delegatee = null;
        pm = null;
    }

    public ProgressReportingDelegatingCodec(final FeatureCodec<A, B> delegatee, final double secondsBetweenUpdates){
        if ( secondsBetweenUpdates <= 0.0 ) {
            throw new IllegalArgumentException("secondsBetweenUpdates must be > 0.0");
        }
        this.delegatee = delegatee;
        this.pm = new ProgressMeter(secondsBetweenUpdates);
    }

    @Override
    public Feature decodeLoc(final B b) throws IOException {
        if (delegatee == null) {
            throw new IllegalStateException("this codec cannot be used without a delegatee.");
        }
        if (!pm.started()) {
            pm.start();
        }
        final Feature f = delegatee.decodeLoc(b);
        pm.update(f);
        return f;
    }

    @Override
    public A decode(final B b) throws IOException {
        if (delegatee == null) {
            throw new IllegalStateException("this codec cannot be used without a delegatee.");
        }
        if (!pm.started()) {
            pm.start();
        }

        final A result = delegatee.decode(b);
        pm.update(result);
        return result;
    }

    @Override
    public FeatureCodecHeader readHeader(final B b) throws IOException {
        if (delegatee == null) {
            throw new IllegalStateException("this codec cannot be used without a delegatee.");
        }
        return delegatee.readHeader(b);
    }

    @Override
    public Class<A> getFeatureType() {
        if (delegatee == null) {
            throw new IllegalStateException("this codec cannot be used without a delegatee.");
        }
        return delegatee.getFeatureType();
    }

    @Override
    public B makeSourceFromStream(final InputStream bufferedInputStream) {
        if (delegatee == null) {
            throw new IllegalStateException("this codec cannot be used without a delegatee.");
        }
        return delegatee.makeSourceFromStream(bufferedInputStream);
    }

    @Override
    public LocationAware makeIndexableSourceFromStream(final InputStream bufferedInputStream) {
        if (delegatee == null) {
            throw new IllegalStateException("this codec cannot be used without a delegatee.");
        }
        return delegatee.makeIndexableSourceFromStream(bufferedInputStream);
    }

    @Override
    public boolean isDone(final B b) {
        if (delegatee == null) {
            throw new IllegalStateException("this codec cannot be used without a delegatee.");
        }
        final boolean done = delegatee.isDone(b);
        if (done){
            pm.stop();
        }
        return done;
    }

    @Override
    public void close(final B b) {
        if (delegatee == null) {
            throw new IllegalStateException("this codec cannot be used without a delegatee.");
        }
        delegatee.close(b);
    }

    @Override
    public boolean canDecode(final String path) {
        //If there's no delegatee then we're going to say no to all questions here
        return delegatee != null && delegatee.canDecode(path);
    }

    public FeatureCodec<A, B> getDelegatee() {
        return delegatee;
    }

    @Override
    public TabixFormat getTabixFormat() {
        if (delegatee == null) {
            throw new IllegalStateException("this codec cannot be used without a delegatee.");
        }
        return delegatee.getTabixFormat();
    }
}
