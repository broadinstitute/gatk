package org.broadinstitute.hellbender.utils.logging;

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
    // if false, then we don't print any message.
    private static boolean enabled = true;
    private Logger optLogger = null;
    private final String id;

    public BunnyLog() {
        this.id = "" + new Random().nextLong();
    }

    public BunnyLog(Logger l) {
        this();
        optLogger = l;
    }

    /**
     * Pass false to disable all logging via bunnylog
     * (it's enabled by default).
     */
    public static void setEnabled(boolean enabled) {
        BunnyLog.enabled = enabled;
    }

    public String start(String name) {
        String ret = bunny + " START " + id + " " + name;
        if (enabled && null!=optLogger) optLogger.info(ret);
        return ret;
    }

    public String stepEnd(String stepName) {
        String ret = bunny + " STEPEND " + id + " " + stepName;
        if (enabled && null!=optLogger) optLogger.info(ret);
        return ret;
    }

    public String end() {
        String ret = bunny + " END " + id;
        if (enabled && null!=optLogger) optLogger.info(ret);
        return ret;
    }
}
