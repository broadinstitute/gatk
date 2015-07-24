package org.broadinstitute.hellbender.dev;

import java.util.Random;
import org.apache.logging.log4j.Logger;

/**
 * The "bunny" log format:
 *
 * =[**]= START <id> <name>
 * =[**]= STEPEND <id> <step_name>
 * =[**]= END <id>
 *
 * The functions here create an id for you and keep track of it, and format the various strings,
 * sending it to a logger if you provided one.
 */
public final class BunnyLog {
    // to mark the log entries in the "bunny" format
    public static final String bunny = "=[**]=";
    private Logger optLogger = null;
    private final String id;

    public BunnyLog() {
        this.id = "" + new Random().nextLong();
    }

    public BunnyLog(Logger l) {
        this();
        optLogger = l;
    }

    public String start(String name) {
        String ret = bunny + " START " + id + " " + name;
        if (null!=optLogger) optLogger.info(ret);
        return ret;
    }

    public String stepEnd(String stepName) {
        String ret = bunny + " STEPEND " + id + " " + stepName;
        if (null!=optLogger) optLogger.info(ret);
        return ret;
    }

    public String end() {
        String ret = bunny + " END " + id;
        if (null!=optLogger) optLogger.info(ret);
        return ret;
    }
}
