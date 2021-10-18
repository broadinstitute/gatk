package org.broadinstitute.hellbender.tools.copynumber.arguments;

import com.google.common.collect.ImmutableList;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.hellbender.tools.copynumber.GermlineCNVCaller;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class GermlineCallingArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final String P_ALT_LONG_NAME = "p-alt";
    public static final String P_ACTIVE_LONG_NAME = "p-active";
    public static final String CNV_COHERENCE_LENGTH_LONG_NAME = "cnv-coherence-length";
    public static final String CLASS_COHERENCE_LENGTH_LONG_NAME = "class-coherence-length";
    public static final String MAX_COPY_NUMBER_LONG_NAME = "max-copy-number";

    // these model parameters will be extracted from provided model in CASE mode
    private static final List<String> HIDDEN_ARGS_CASE_MODE = ImmutableList.of(
            P_ACTIVE_LONG_NAME,
            CLASS_COHERENCE_LENGTH_LONG_NAME);

    @Argument(
            doc = "Total prior probability of alternative copy-number states (the reference copy-number " +
                    "is set to the contig integer ploidy)",
            fullName = P_ALT_LONG_NAME,
            minValue = 0.,
            maxValue = 1.,
            optional = true
    )
    private double pAlt = 1e-6;

    @Argument(
            doc = "Prior probability of treating an interval as CNV-active (in a CNV-active domains, all " +
                    "copy-number states are equally likely to be called).",
            fullName = P_ACTIVE_LONG_NAME,
            minValue = 0.,
            maxValue = 1.,
            optional = true
    )
    private double pActive = 1e-2;

    @Argument(
            doc = "Coherence length of CNV events (in the units of bp).",
            fullName = CNV_COHERENCE_LENGTH_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double cnvCoherenceLength = 10000.0;

    @Argument(
            doc = "Coherence length of CNV-active and CNV-silent domains (in the units of bp).",
            fullName = CLASS_COHERENCE_LENGTH_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double classCoherenceLength = 10000.0;

    @Argument(
            doc = "Highest allowed copy-number state.",
            fullName = MAX_COPY_NUMBER_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int maxCopyNumber = 5;

    public List<String> generatePythonArguments(final GermlineCNVCaller.RunMode runMode) {
        final List<String> arguments = new ArrayList<>(Arrays.asList(
                String.format("--p_alt=%e", pAlt),
                String.format("--cnv_coherence_length=%e", cnvCoherenceLength),
                String.format("--max_copy_number=%d", maxCopyNumber)));
        if (runMode == GermlineCNVCaller.RunMode.COHORT) {
            arguments.addAll(Arrays.asList(
                    String.format("--p_active=%f", pActive),
                    String.format("--class_coherence_length=%f", classCoherenceLength)));
        }
        return arguments;
    }

    public void validate(final CommandLineArgumentParser clpParser, final GermlineCNVCaller.RunMode runMode) {
        if (runMode == GermlineCNVCaller.RunMode.CASE)
            HIDDEN_ARGS_CASE_MODE.forEach(a -> Utils.validateArg(
                    !clpParser.getNamedArgumentDefinitionByAlias(a).getHasBeenSet(),
                    String.format("Argument '--%s' cannot be set in the CASE mode.", a)));
        ParamUtils.isPositive(cnvCoherenceLength,
                String.format("Coherence length of CNV events (%s) must be positive.",
                        CNV_COHERENCE_LENGTH_LONG_NAME));
        ParamUtils.isPositive(classCoherenceLength,
                String.format("Coherence length of CNV class domains (%s) must be positive.",
                        CLASS_COHERENCE_LENGTH_LONG_NAME));
        ParamUtils.isPositive(maxCopyNumber,
                String.format("Highest allowed copy-number (%s) must be positive.",
                        MAX_COPY_NUMBER_LONG_NAME));
    }
}
