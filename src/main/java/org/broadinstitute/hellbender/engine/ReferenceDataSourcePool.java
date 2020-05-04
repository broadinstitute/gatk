package org.broadinstitute.hellbender.engine;

import org.apache.commons.pool.PoolableObjectFactory;
import org.apache.commons.pool.impl.GenericObjectPool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.AutoCloseableReference;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

class ReferenceDataSourcePool extends DataSourcePool<ReferenceDataSource> {

    ReferenceDataSourcePool(final Path fastaPath) {
        super(new Factory(fastaPath));
    }

    private static class Factory implements PoolableObjectFactory<ReferenceDataSource> {

        private final Path path;

        private Factory(final Path fastaPath) {
            this.path = fastaPath;
        }

        @Override
        public ReferenceDataSource makeObject() throws Exception {
            return ReferenceDataSource.of(path);
        }

        @Override
        public void destroyObject(ReferenceDataSource dataSource) {
            dataSource.close();
        }

        @Override
        public boolean validateObject(ReferenceDataSource dataSource) {
            return true;
        }

        @Override
        public void activateObject(ReferenceDataSource dataSource) {
            // do nothing.
        }

        @Override
        public void passivateObject(ReferenceDataSource dataSource) {
            // do nothing.
        }
    }
}
