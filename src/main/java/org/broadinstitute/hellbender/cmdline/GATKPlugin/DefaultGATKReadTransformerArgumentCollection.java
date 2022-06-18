package org.broadinstitute.hellbender.cmdline.GATKPlugin;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.ReadTransformerArgumentDefinitions;

import java.util.ArrayList;
import java.util.List;

/**
 * Default {@link GATKReadTransformerArgumentCollection} applied in GATK for optional read transformers in the command line.
 * It contains arguments that allow the user to:
 *
 * - Provide a list of read transformers to apply.
 * - Disable some and/or all read transformers.
 *
 * @author Dror Kessler (dror27)
 */
public class DefaultGATKReadTransformerArgumentCollection extends GATKReadTransformerArgumentCollection {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = ReadTransformerArgumentDefinitions.READ_TRANSFORMER_LONG_NAME,
            shortName = ReadTransformerArgumentDefinitions.READ_TRANSFORMER_SHORT_NAME,
            doc="Read transformers to be applied before analysis", optional=true, common = true)
    public final List<String> userEnabledReadTransformerNames = new ArrayList<>(); // preserve order

    @Argument(fullName = ReadTransformerArgumentDefinitions.DISABLE_READ_TRANSFORMER_LONG_NAME,
            shortName = ReadTransformerArgumentDefinitions.DISABLE_READ_TRANSFORMER_SHORT_NAME,
            doc="Read transformers to be disabled before analysis", optional=true, common = true)
    public final List<String> userDisabledReadTransformerNames = new ArrayList<>();

    @Advanced
    @Argument(fullName = ReadTransformerArgumentDefinitions.DISABLE_TOOL_DEFAULT_READ_TRANSFORMERS,
            shortName = ReadTransformerArgumentDefinitions.DISABLE_TOOL_DEFAULT_READ_TRANSFORMERS,
            doc = "Disable all tool default read transformers (WARNING: many tools will not function correctly without their default read transformers on)", common = true, optional = true)
    public boolean disableToolDefaultReadTransformers = false;

    /** Returns the list with the read transformers provided by the user, preserving the order. */
    @Override
    public List<String> getUserEnabledReadTransformerNames() {
        return userEnabledReadTransformerNames;
    }

    /** Returns the set of transformers disabled by the user. */
    @Override
    public List<String> getUserDisabledReadTransformerNames() {
        return userDisabledReadTransformerNames;
    }

    /** {@inheritDoc}. */
    @Override
    public boolean getDisableToolDefaultReadTransformers() {
        return disableToolDefaultReadTransformers;
    }

}
