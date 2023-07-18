package org.broadinstitute.hellbender.utils.help;

import htsjdk.samtools.metrics.MetricBase;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DefaultDocWorkUnitHandler;
import org.broadinstitute.barclay.help.DocWorkUnit;
import org.broadinstitute.barclay.help.HelpDoclet;
import org.broadinstitute.barclay.help.scanners.JavaLanguageModelScanners;

import javax.lang.model.element.Element;
import javax.lang.model.element.ElementKind;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
    private final static String PICARD_METRICS_TEMPLATE_NAME = "metrics.template.html";

    private final static String WORK_UNIT_SUMMARY_KEY = "summary";
    private final static String METRICS_MAP_ENTRY_KEY = "metrics";
    private final static String METRICS_MAP_NAME_KEY = "name";
    private final static String METRICS_MAP_SUMMARY_KEY = "summary";

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
     * @return the name of the freemarker template to be used for the class being documented.
     * Must reside in the folder passed to the Barclay Doclet via the "-settings-dir" parameter to
     * Javadoc.
     */
    @Override
    public String getTemplateName(final DocWorkUnit workUnit) {
        Class<?> clazz = workUnit.getClazz();
        if (MetricBase.class.isAssignableFrom(clazz)) {
            return PICARD_METRICS_TEMPLATE_NAME;
        } else {
            return GATK_FREEMARKER_TEMPLATE_NAME;
        }
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
        } else if (MetricBase.class.isAssignableFrom(toolClass)) {
            currentWorkUnit.setProperty(WORK_UNIT_SUMMARY_KEY, currentWorkUnit.getSummary());
            final List<Map<String, String>> workUnitMetricsList = new ArrayList<>();
            currentWorkUnit.setProperty(METRICS_MAP_ENTRY_KEY, workUnitMetricsList);
            final Field[] fields = currentWorkUnit.getClazz().getFields();
            for (final Field field : fields) {
                if (Modifier.isPublic(field.getModifiers())) {
                    final Element fieldElement = JavaLanguageModelScanners.getElementForField(
                            getDoclet().getDocletEnv(),
                            currentWorkUnit.getDocElement(),
                            field,
                            ElementKind.FIELD
                    );
                    if (fieldElement != null) {
                        final String docComment = JavaLanguageModelScanners.getDocComment(getDoclet().getDocletEnv(), fieldElement);
                        final Map<String, String> metricsFields = new HashMap<>();
                        metricsFields.put(METRICS_MAP_NAME_KEY, field.getName());
                        metricsFields.put(METRICS_MAP_SUMMARY_KEY, docComment);
                        workUnitMetricsList.add(metricsFields);
                    }
                }
            }
        }
    }

}
