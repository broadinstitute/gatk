package org.broadinstitute.hellbender.engine.spark;


import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Log;
import org.apache.logging.log4j.Level;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Command line arguments needed for configuring a spark context
 */
public final class SparkCommandLineArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final String SPARK_MASTER_LONG_NAME = "spark-master";
    public static final String SPARK_VERBOSITY_LONG_NAME = "spark-verbosity";

    @Argument(
            doc="URL of the Spark Master to submit jobs to when using the Spark pipeline runner.",
            fullName = SPARK_MASTER_LONG_NAME,
            optional = true)
    private String sparkMaster = SparkContextFactory.DEFAULT_SPARK_MASTER;

    @Argument(
            doc = "Spark properties to set on the Spark context in the format <property>=<value>",
            fullName = StandardArgumentDefinitions.SPARK_PROPERTY_NAME,
            optional = true
    )
    final List<String> sparkProperties = new ArrayList<>();

    @Argument(
            doc="Spark verbosity. Overrides --" + StandardArgumentDefinitions.VERBOSITY_NAME + " for Spark-generated logs only. Possible values: {ALL, DEBUG, INFO, WARN, ERROR, FATAL, OFF, TRACE}",
            fullName = SPARK_VERBOSITY_LONG_NAME,
            optional = true)
    private String sparkVerbosity = null;

    public Map<String,String> getSparkProperties(){
        final Map<String, String> propertyMap = new LinkedHashMap<>();
        for( String property: sparkProperties) {
            final String[] splits = property.split("=");
            if (splits.length != 2 || splits[0].isEmpty() || splits[1].isEmpty()) {
                throw new CommandLineException.BadArgumentValue(StandardArgumentDefinitions.SPARK_PROPERTY_NAME, property, "Expected a value of the form spark.property.name=value");
            } else {
                propertyMap.put(splits[0], splits[1]);
            }
        }
        return propertyMap;
    }

    public String getSparkMaster() {
        return sparkMaster;
    }

    /**
     * Returns the Spark log level for the argument set. This is simply sparkVerbosity
     * if it was specified. Otherwise, it returns the log level corresponding to the
     * provided tool (htsjdk) verbosity.
     *
     * @param  toolVerbosity  Current tool's htsjdk log level
     * @return      Spark log level String
     */
    public String getSparkVerbosity(final Log.LogLevel toolVerbosity) {
        Utils.nonNull(toolVerbosity, "Tool verbosity cannot be null");
        if (sparkVerbosity != null) return sparkVerbosity;
        if (toolVerbosity.equals(Log.LogLevel.DEBUG)) {
            return Level.DEBUG.name();
        }
        if (toolVerbosity.equals(Log.LogLevel.INFO)) {
            return Level.INFO.name();
        }
        if (toolVerbosity.equals(Log.LogLevel.WARNING)) {
            return Level.WARN.name();
        }
        if (toolVerbosity.equals(Log.LogLevel.ERROR)) {
            return Level.ERROR.name();
        }
        throw new IllegalStateException("Unknown tool verbosity: " + toolVerbosity.name());
    }

    @VisibleForTesting
    public void setSparkVerbosity(String level) {
        sparkVerbosity = level;
    }
}
