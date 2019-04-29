package org.broadinstitute.hellbender.tools.walkers.varianteval;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.AlleleFrequency.StratifyingScale;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.io.Resource;

import java.util.Arrays;
import java.util.Collections;


/**
 * Created by skwalker on 4/19/19.
 */
public class AlleleFrequencyQC extends VariantEval {

    private String R_SCRIPT = "plotAlleleFrequencyQC.R";

    @Argument(shortName = "sample-name",
            doc="Sample name to be used in metric file")
    protected String sample_name; // can we get this from the vcf ?


    @Override
    public void onTraversalStart() {
        NO_STANDARD_MODULES = true;
        MODULES_TO_USE = Collections.singletonList("PVariantEvaluator");
        keepSitesWithAC0 = true;
        NO_STANDARD_STRATIFICATIONS = true;
        STRATIFICATIONS_TO_USE = Arrays.asList("AlleleFrequency", "Filter");
        AFScale = StratifyingScale.LOGARITHMIC;
        useCompAFStratifier = true;

        super.onTraversalStart();
    }

    @Override
    public Object  onTraversalSuccess() {

        super.onTraversalSuccess();

        // need the file returned from variant eval in order to run the plotting stuff
        final RScriptExecutor executer = new RScriptExecutor();
        executer.addScript(new Resource(R_SCRIPT, AlleleFrequencyQC.class));
        executer.addArgs(outFile.getAbsoluteFile(), outFile.getAbsolutePath() + sample_name, sample_name);
        executer.exec();

        return null;
    }



}
