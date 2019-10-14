package org.broadinstitute.hellbender.utils.help;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;

import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocWorkUnit;
import org.broadinstitute.barclay.help.HelpDoclet;
import org.broadinstitute.barclay.help.WDLWorkUnitHandler;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKPathSpecifier;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.*;

// Note: WDL Gen doesn't handle arguments that accept tagged argument values

/**
 * The GATK WDL work unit handler. Its main task is to convert the types for all arguments for a given work
 * unit (tool) from Java types to WDL-compatible types by updating the freemarker map with the transformed types.
 *
 * NOTE: Methods in this class are intended to be called by Gradle/Javadoc only, and should not be called
 * by methods that are used by the GATK runtime, as this class assumes a dependency on com.sun.javadoc classes
 * which may not be present.
 */
public class GATKWDLWorkUnitHandler extends WDLWorkUnitHandler {

    private final static String GATK_FREEMARKER_TEMPLATE_NAME = "wdlToolTemplate.wdl.ftl";

    // Map of Java argument types that the WDL generator knows how to convert to a WDL type, along with the
    // corresponding string substitution that needs to be run on the (Barclay-generated) string that describes
    // the type. From a purely string perspective, some of these transforms are no-ops in that no actual
    // conversion is required because the type names are identical in Java and WDL (i.e, File->File or
    // String->String), but they're included here for completeness, and to document the allowed type transitions.
    private final static Map<Class<?>, ImmutablePair<String, String>> javaToWDLTypeMap =
            new HashMap<Class<?>, ImmutablePair<String, String>>() {
                private static final long serialVersionUID = 1L;
                {
                    // GATK-specific File Types
                    put(GATKPathSpecifier.class, new ImmutablePair<>(GATKPathSpecifier.class.getSimpleName(), "File"));
                    // FeatureInputs require special handling to account for the generic type param(s)
                    put(FeatureInput.class, new ImmutablePair<>(FeatureInput.class.getSimpleName(), "File"));
            }
        };

    public GATKWDLWorkUnitHandler(final HelpDoclet doclet) {
        super(doclet);
    }

    /**
     * @param workUnit the DocWorkUnit object being processed
     * @return the name of a the freemarker template to be used for the class being documented.
     * Must reside in the folder passed to the Barclay Doclet via the "-settings-dir" parameter to
     * Javadoc.
     */
    @Override
    public String getTemplateName(final DocWorkUnit workUnit) { return GATK_FREEMARKER_TEMPLATE_NAME; }

    /**
     * Return the flat filename (no paths) that the handler would like the Doclet to
     * write out the documentation for workUnit.
     * @param workUnit
     * @return the name of the destination file to which documentation output will be written
     */
    @Override
    public String getDestinationFilename(final DocWorkUnit workUnit) {
        return workUnit.getClazz().getSimpleName() + ".wdl";
    }

    /**
     * Returns the JSON output file name.
     */
    @Override
    public String getJSONFilename(final DocWorkUnit workUnit) {
        return workUnit.getClazz().getSimpleName() + "Inputs.json";
    }

    /**
     * Given a Java class representing the underlying field  type of an argument, and a human readable doc type,
     * return a String with the corresponding WDL type.
     *
     * @param argumentClass the Class for the underlying field of the argument being converted
     * @param docType a string representing the human readable type assigned by the Barclay doc system
     * @param contextMessage a message describing the context for this argument, used in error reporting
     * @return the docType string transformed to the corresponding WDL type
     */
    @Override
    protected String convertJavaTypeToWDLType(final Class<?> argumentClass, final String docType, final String contextMessage) {
        String convertedWDLType;
        if (FeatureInput.class.isAssignableFrom(argumentClass)) {
            if (!docType.contains(FeatureInput.class.getSimpleName())) {
                throw new GATKException(
                        String.format(
                                "Don't know how to generate a WDL type for %s in work unit %s",
                                argumentClass,
                                contextMessage));
            }
            final Pair<String, String> typeConversionPair = getJavaWDLTypeConversion(argumentClass);
            convertedWDLType = docType.replaceFirst("FeatureInput\\[[a-zA-Z0-9?]+\\]", typeConversionPair.getValue());
        } else {
            return super.convertJavaTypeToWDLType(argumentClass, docType, contextMessage);
        }
        return convertedWDLType;
    }

    /**
     * Given an argument class, return a String pair representing the string that should be replaced (the Java type),
     * and the string to substitute (the corresponding WDL type), i.e., for an argument with type Java Integer.class,
     * return the Pair ("Integer", "Int") to convert from the Java type to the corresponding WDL type.
     * @param argumentClass Class of the argument being converter
     * @return a String pair representing the original and replacement type text, or null if no conversion is available
     */
    @Override
    protected Pair<String, String> getJavaWDLTypeConversion(final Class<?> argumentClass) {
        Pair<String, String> conversion = javaToWDLTypeMap.get(argumentClass);
        return conversion == null ? super.getJavaWDLTypeConversion(argumentClass) : conversion;
    }

    /**
     * Add any custom freemarker bindings discovered via custom javadoc tags. Subclasses can override this to
     * provide additional custom bindings.
     *
     * @param currentWorkUnit the work unit for the feature being documented
     */
    @Override
    protected void addCustomBindings(final DocWorkUnit currentWorkUnit) {
        super.addCustomBindings(currentWorkUnit);

        // Picard tools use the summary line for the long overview section, so extract that
        // from Picard tools only, and put it in the freemarker map.
        Class<?> toolClass = currentWorkUnit.getClazz();
        if (picard.cmdline.CommandLineProgram.class.isAssignableFrom(toolClass)) {
            final CommandLineProgramProperties clpProperties = currentWorkUnit.getCommandLineProperties();
            currentWorkUnit.setProperty("picardsummary", clpProperties.summary());
        }
    }

}
