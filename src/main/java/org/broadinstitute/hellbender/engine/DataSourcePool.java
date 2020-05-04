package org.broadinstitute.hellbender.engine;

import org.apache.commons.pool.PoolableObjectFactory;
import org.apache.commons.pool.impl.GenericObjectPool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.AutoCloseableReference;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

class DataSourcePool<T> extends GenericObjectPool<T> implements AutoCloseable {

    DataSourcePool(final PoolableObjectFactory<T> factory) {
        super(factory);
        setWhenExhaustedAction(WHEN_EXHAUSTED_GROW);
    }

    public AutoReturn borrowAutoReturn() {
        return new AutoReturn(borrowObject());
    }

    @Override
    public T borrowObject() {
        try {
            return super.borrowObject();
        } catch (Exception e) {
            throw new GATKException(e.getMessage(), e);
        }
    }

    @Override
    public void returnObject(T dataSource) {
        try {
            super.returnObject(dataSource);
        } catch (Exception e) {
            throw new GATKException(e.getMessage(), e);
        }
    }

    public final class AutoReturn extends AutoCloseableReference<T> {

        private AutoReturn(final T subject) {
            super(subject);
        }


        @Override
        protected void close(T subject) {
            returnObject(subject);
        }
    }

    public void close() {
        try {
            super.close();
        } catch (final Exception ex) {
            throw new GATKException("exception enoutred", ex);
        }
    }

}
