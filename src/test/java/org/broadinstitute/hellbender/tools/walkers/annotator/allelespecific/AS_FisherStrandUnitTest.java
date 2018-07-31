package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;

import static org.testng.Assert.*;

/**
 * This class exists in order to run GATK3 concordance tests using the framework of ReducibleAnnotationBaseTest
 * Created by emeryj on 8/11/17.
 */
public class AS_FisherStrandUnitTest extends ReducibleAnnotationBaseTest {

    @Override
    protected List<Annotation> getAnnotationsToUse() {
        return Collections.singletonList(new AS_FisherStrand());
    }

    @Override
    protected String getRawKey() {
        return GATKVCFConstants.AS_FISHER_STRAND_KEY;
    }

    @Override
    protected String getKey() {
        return GATKVCFConstants.AS_FISHER_STRAND_KEY;
    }
}