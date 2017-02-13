package org.broadinstitute.hellbender.utils.help;

import org.broadinstitute.barclay.help.DefaultDocWorkUnitHandler;
import org.broadinstitute.barclay.help.DocWorkUnit;

import org.broadinstitute.barclay.help.HelpDoclet;

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
}
