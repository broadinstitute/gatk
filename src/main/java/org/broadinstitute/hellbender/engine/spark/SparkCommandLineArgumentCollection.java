package org.broadinstitute.hellbender.engine.spark;


import org.apache.logging.log4j.Level;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

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
            shortName = StandardArgumentDefinitions.SPARK_PROPERTY_NAME,
            optional = true
    )
    final List<String> sparkProperties = new ArrayList<>();

    @Argument(
            doc="Spark verbosity (ALL, DEBUG, INFO, WARN, ERROR, FATAL, OFF, TRACE)",
            fullName = SPARK_VERBOSITY_LONG_NAME,
            optional = true)
    private String sparkVerbosity = Level.INFO.name();

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
    public String getSparkVerbosity() {
        return sparkVerbosity;
    }
}
