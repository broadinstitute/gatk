package org.broadinstitute.hellbender.engine;

import org.apache.commons.pool.PoolableObjectFactory;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

class ReadsDataSourcePool extends DataSourcePool<ReadsDataSource> {

    ReadsDataSourcePool(final List<Path> readPaths) {
        super(new Factory(readPaths));
        setWhenExhaustedAction(WHEN_EXHAUSTED_GROW);
    }


    private static class Factory implements PoolableObjectFactory<ReadsDataSource> {

        private final List<Path> paths;

        private Factory(final List<Path> paths) {
            this.paths = new ArrayList<>(paths);
        }

        @Override
        public ReadsDataSource makeObject() {
            return new ReadsDataSource(paths);
        }

        @Override
        public void destroyObject(ReadsDataSource dataSource) {
            dataSource.close();
        }

        @Override
        public boolean validateObject(ReadsDataSource dataSource) {
            return true;
        }

        @Override
        public void activateObject(ReadsDataSource dataSource) {
            // do nothing.
        }

        @Override
        public void passivateObject(ReadsDataSource dataSource) {
            // do nothing.
        }
    }
}
