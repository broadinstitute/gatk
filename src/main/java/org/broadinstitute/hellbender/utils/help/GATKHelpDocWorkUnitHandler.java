package org.broadinstitute.hellbender.utils.help;

import org.broadinstitute.barclay.argparser.CommandLinePluginDescriptor;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.NamedArgumentDefinition;
import org.broadinstitute.barclay.help.DefaultDocWorkUnitHandler;
import org.broadinstitute.barclay.help.DocWorkUnit;

import org.broadinstitute.barclay.help.HelpDoclet;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKAnnotationPluginDescriptor;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import java.util.*;
import java.util.stream.Collectors;

/**
 * The GATK Documentation work unit handler class that is the companion to GATKHelpDoclet.
 *
 * NOTE: Methods in this class are intended to be called by Gradle/Javadoc only, and should not be called
 * by methods that are used by the GATK runtime, as this class assumes a dependency on com.sun.javadoc classes
 * which may not be present.
 */
public class GATKHelpDocWorkUnitHandler extends DefaultDocWorkUnitHandler {

    private final static String GATK_JAVADOC_TAG_PREFIX = "GATK"; // prefix for custom javadoc tags used by GATK

    private final static String GATK_FREEMARKER_TEMPLATE_NAME = "generic.template.html";

    public GATKHelpDocWorkUnitHandler(final HelpDoclet doclet) {
        super(doclet);
    }
    /**
     * @return Prefix for custom GATK tags that should be lifted from the javadoc and stored in the
     * FreeMarker map. These will be available in the template returned by {@link #getTemplateName}.
     */
    @Override
    protected String getTagFilterPrefix() { return GATK_JAVADOC_TAG_PREFIX; }

    /**
     * @param workUnit the classdoc object being processed
     * @return the name of a the freemarker template to be used for the class being documented.
     * Must reside in the folder passed to the Barclay Doclet via the "-settings-dir" parameter to
     * Javadoc.
     */
    @Override
    public String getTemplateName(final DocWorkUnit workUnit) { return GATK_FREEMARKER_TEMPLATE_NAME; }


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

    @Override
    protected void addDefaultPlugins(
            final DocWorkUnit currentWorkUnit,
            final List<? extends CommandLinePluginDescriptor<?>> pluginDescriptors) {
        super.addDefaultPlugins(currentWorkUnit, pluginDescriptors);

        final String readFilterDescriptorName = new GATKReadFilterPluginDescriptor(
                Collections.emptyList()).getDisplayName();
        final String annotationDescriptorName = new GATKAnnotationPluginDescriptor(
                Collections.emptyList(),
                Collections.emptyList()).getDisplayName();

        // for the --read-filters and --annotations arguments, we need to artificially present the default
        // plugin values that are programmatically set by tools as if they were "default" values for the corresponding
        // arguments so they show up in the doc for each tool
        pluginDescriptors.forEach(
                descriptor -> {
                    if (descriptor.getDisplayName().equals(readFilterDescriptorName)) {
                        propagatePluginDefaults(
                                currentWorkUnit,
                                readFilterDescriptorName,
                                ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME);
                    } else if (descriptor.getDisplayName().equals(annotationDescriptorName)) {
                        propagatePluginDefaults(
                                currentWorkUnit,
                                annotationDescriptorName,
                                StandardArgumentDefinitions.ANNOTATION_LONG_NAME);
                    }
                });

    }

    // add the default instances for plugin "descriptorName" as default values for "targetArgumentName"
    @SuppressWarnings("unchecked")
    private void propagatePluginDefaults(
            final DocWorkUnit currentWorkUnit,
            final String descriptorName,
            final String targetArgumentName) {
        final List<String> defaultReadFilterNames = new ArrayList<>();

        final HashSet<HashMap<String, Object>> defaultsForPlugins =
                (HashSet<HashMap<String, Object>>) currentWorkUnit.getProperty(descriptorName);
        if (defaultsForPlugins != null) {
            for (final HashMap<String, Object> pluginMap : defaultsForPlugins) {
                defaultReadFilterNames.add((String) pluginMap.get("name"));
            }
        }
        final Map<String, List<Map<String, Object>>> argsMap =
                (Map<String, List<Map<String, Object>>>) currentWorkUnit.getProperty("arguments");
        final List<Map<String, Object>> readFilterArgList =
                argsMap.get("all").stream().filter(
                        m -> m.containsKey("name") && m.get("name").equals("--" + targetArgumentName)
                ).collect(Collectors.toList());
        if (readFilterArgList.size() != 1) {
            throw new IllegalStateException(String.format("Can't find argument %s for descriptor %s in %s",
                    targetArgumentName,
                    descriptorName,
                    currentWorkUnit.getClazz()));
        }
        final Map<String, Object> readFilterArg = readFilterArgList.get(0);
        final String oldDefaultValue = (String) readFilterArg.get("defaultValue");
        if (oldDefaultValue != null &&
                oldDefaultValue.length() != 0 && !oldDefaultValue.equals("") && !oldDefaultValue.equals("[]")) {
            throw new IllegalStateException(String.format("%s argument property for argument %s for descriptor %s in %s is already populated",
                    "defaultValue",
                    targetArgumentName,
                    currentWorkUnit.getClazz()));
        }
        readFilterArg.put("defaultValue", defaultReadFilterNames.stream().collect(Collectors.joining(", ")));
    }

}
