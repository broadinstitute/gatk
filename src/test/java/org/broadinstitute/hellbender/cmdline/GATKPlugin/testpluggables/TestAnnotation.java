package org.broadinstitute.hellbender.cmdline.GATKPlugin.testpluggables;

import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;

/**
 * This is a test annotation that is used to test if the gatk config file properly controls annotation loading
 * see {@link org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKAnnotationPluginDescriptorUnitTest}
 *
 * This should not be discoverable normally since it's not in the expected annotation package.
 */
public class TestAnnotation implements Annotation {
}
