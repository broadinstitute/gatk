package org.broadinstitute.hellbender.utils.logging;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.simple.SimpleLogger;
import org.apache.logging.log4j.util.PropertiesUtil;
import org.testng.Assert;
import org.testng.annotations.Test;
import java.util.Properties;


public class OneShotLoggerTest {

    private class StorageLogger extends SimpleLogger {
        private static final long serialVersionUID = 1L;
        int timesLogged = 0;

        public StorageLogger() {
            super("", Level.ALL, false, false, false, false, "", null, new PropertiesUtil(new Properties()), null);
        }

        @Override
        public void warn(String message) {
            timesLogged++;
        }
    }

    @Test
    public void testStart() throws Exception {
        OneShotLogger log = new OneShotLogger(this.getClass());
        log.logger = new StorageLogger();
        log.warn("Foo");
        Assert.assertEquals(((StorageLogger)log.logger).timesLogged , 1);
        log.warn("Bar");
        log.warn("Baz");
        Assert.assertEquals(((StorageLogger)log.logger).timesLogged , 1);
    }
}