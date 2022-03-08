package org.broadinstitute.hellbender.tools.cfCNV;

import breeze.stats.distributions.NegativeBinomial;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;

import java.io.File;
import java.util.List;

@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = CopyNumberProgramGroup.class // sato: change
)
@DocumentedFeature
public class CfCNVDetect extends GATKTool {
    // sato: define name variables in CopyNumberStandardArgument
    @Argument(doc = "", fullName = "", minElements = 1)
    private File inputDenoisedCopyRatiosFile;

    @Argument(doc = "",fullName = "", minElements = 1)
    private File inputAllelicCountsFile;

    @Override
    public void traverse() {
        final CopyRatioCollection crc = new CopyRatioCollection(inputDenoisedCopyRatiosFile);
        final AllelicCountCollection acc = new AllelicCountCollection(inputAllelicCountsFile);

        // Identify copy neutral spots

        // Fit negative binomial....or Gaussian?

        // [ Comment: need to get used to using interval intersection in python ]

        // Want to be able to plot stuff...hmm, ok, back to python.
        final NegativeBinomial nb = new NegativeBinomial(1.0, 1.0, null);
    }
}
