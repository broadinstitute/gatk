package org.broadinstitute.hellbender.engine.spark;


import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Command line arguments needed for configuring a spark context
 */
public final class SparkCommandLineArgumentCollection implements ArgumentCollectionDefinition {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "whether to use the Spark implementation (default is not to use Spark)", shortName = "spark", fullName = "spark", optional = true)
    public boolean useSpark = false;

    @Argument(fullName = "sparkMaster", doc="URL of the Spark Master to submit jobs to when using the Spark pipeline runner.", optional = true)
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
            final String[] splits = property.split("=");
            if (splits.length != 2 || splits[0].isEmpty() || splits[1].isEmpty()) {
                throw new UserException.BadArgumentValue(StandardArgumentDefinitions.SPARK_PROPERTY_NAME, property, "Expected a value of the form spark.property.name=value");
            } else {
                propertyMap.put(splits[0], splits[1]);
            }
        }
        return propertyMap;
    }

    public String getSparkMaster() {
        return sparkMaster;
    }

}
