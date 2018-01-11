package org.broadinstitute.hellbender.engine.spark;


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

    @Argument(fullName = "spark-master", doc="URL of the Spark Master to submit jobs to when using the Spark pipeline runner.", optional = true)
    private String sparkMaster = SparkContextFactory.DEFAULT_SPARK_MASTER;

    @Argument(
            doc = "spark properties to set on the spark context in the format <property>=<value>",
            shortName = StandardArgumentDefinitions.SPARK_PROPERTY_NAME,
            fullName = StandardArgumentDefinitions.SPARK_PROPERTY_NAME,
            optional = true
    )
    final List<String> sparkProperties = new ArrayList<>();

    public Map<String,String> getSparkProperties(){
        final Map<String, String> propertyMap = new LinkedHashMap<>();
        for( String property: sparkProperties) {
            final List<String> splits = Utils.split(property, '=');
            if (splits.size() != 2 || splits.get(0).isEmpty() || splits.get(1).isEmpty()) {
                throw new CommandLineException.BadArgumentValue(StandardArgumentDefinitions.SPARK_PROPERTY_NAME, property, "Expected a value of the form spark.property.name=value");
            } else {
                propertyMap.put(splits.get(0), splits.get(1));
            }
        }
        return propertyMap;
    }

    public String getSparkMaster() {
        return sparkMaster;
    }

}
