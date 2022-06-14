package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SamReaderFactory;
import org.apache.commons.pool.BasePoolableObjectFactory;
import org.apache.commons.pool.impl.GenericObjectPool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.AutoCloseableReference;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * Pool of {@link ReadsDataSource} instances.
 */
public final class ReadsDataSourcePool extends GenericObjectPool<ReadsPathDataSource> implements AutoCloseable {

    public ReadsDataSourcePool(final List<Path> readPaths, final Path referencePath) {
        super(new Factory(readPaths, referencePath));
        setWhenExhaustedAction(WHEN_EXHAUSTED_GROW);
    }

    /**
     * Returns a reads-data-source wrapped into a reference that when close returns
     * the source back to the pool.
     * @return never {@code null}.
     */
    public AutoCloseableReference<ReadsPathDataSource> borrowAutoReturn() {
        return AutoCloseableReference.of(borrowObject(), this::returnObject);
    }

    @Override
    public ReadsPathDataSource borrowObject() {
        try {
            return super.borrowObject();
        } catch (Exception e) {
            throw new GATKException(e.getMessage(), e);
        }
    }

    @Override
    public void returnObject(final ReadsPathDataSource dataSource) {
        try {
            super.returnObject(dataSource);
        } catch (Exception e) {
            throw new GATKException(e.getMessage(), e);
        }
    }

    public void close() {
        try {
            super.close();
        } catch (final Exception ex) {
            throw new GATKException("exception when closing the pool", ex);
        }
    }

    private static class Factory extends BasePoolableObjectFactory<ReadsPathDataSource> {

        private final List<Path> paths;
        private final Path referencePath;

        private Factory(final List<Path> paths, final Path referencePath) {
            this.paths = paths;
            this.referencePath = referencePath;
        }

        @Override
        public ReadsPathDataSource makeObject() {
            SamReaderFactory factory = SamReaderFactory.makeDefault();
            if (referencePath != null) {
                factory.referenceSequence(referencePath);
            }
            return new ReadsPathDataSource(paths, factory);
        }

        @Override
        public void destroyObject(final ReadsPathDataSource dataSource) {
            dataSource.close();
        }

        @Override
        public boolean validateObject(final ReadsPathDataSource dataSource) {
            return !dataSource.isClosed();
        }
    }
}
