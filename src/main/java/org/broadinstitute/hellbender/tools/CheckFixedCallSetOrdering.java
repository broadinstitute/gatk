package org.broadinstitute.hellbender.tools;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;

@CommandLineProgramProperties(summary = "checkthing", oneLineSummary = "checkthing", programGroup = TestProgramGroup.class)
public class CheckFixedCallSetOrdering extends VariantWalker {
    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        for( String name : variant.getSampleNames()){
            final Genotype genotype = variant.getGenotype(name);
            final String originalSampleName = (String)genotype.getAnyAttribute("SN");
            if(!originalSampleName.equals(name)){
                System.out.println("name:" + name + " != " + originalSampleName);
            }
        }
    }
}
