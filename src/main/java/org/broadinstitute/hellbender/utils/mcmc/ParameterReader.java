package org.broadinstitute.hellbender.utils.mcmc;

import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.File;
import java.io.IOException;
import java.util.AbstractMap;
import java.util.Arrays;
import java.util.Map;

/**
 * Table reader for reading a {@link PosteriorSummary} for each global parameter in a {@link ParameterEnum}.
 *
 * @author Samuel Lee &lt;valentin@broadinstitute.org&gt;
 */
public final class ParameterReader<T extends Enum<T> & ParameterEnum> extends TableReader<Map.Entry<T, PosteriorSummary>> {

    private final Class<T> parameterClass;

    public ParameterReader(final File file, final Class<T> parameterClass) throws IOException {
        super(file);
        this.parameterClass = parameterClass;
    }

    @Override
    protected Map.Entry<T, PosteriorSummary> createRecord(final DataLine dataLine) {
        final String parameterName = dataLine.get(ParameterTableColumn.PARAMETER_NAME);
        final T parameter = Enum.valueOf(parameterClass, parameterName);
        final double center = dataLine.getDouble(ParameterTableColumn.POSTERIOR_MODE);
        final double lower = dataLine.getDouble(ParameterTableColumn.POSTERIOR_LOWER);
        final double upper = dataLine.getDouble(ParameterTableColumn.POSTERIOR_UPPER);
        final DecileCollection deciles = new DecileCollection(Arrays.asList(
                dataLine.getDouble(ParameterTableColumn.POSTERIOR_10),
                dataLine.getDouble(ParameterTableColumn.POSTERIOR_20),
                dataLine.getDouble(ParameterTableColumn.POSTERIOR_30),
                dataLine.getDouble(ParameterTableColumn.POSTERIOR_40),
                dataLine.getDouble(ParameterTableColumn.POSTERIOR_50),
                dataLine.getDouble(ParameterTableColumn.POSTERIOR_60),
                dataLine.getDouble(ParameterTableColumn.POSTERIOR_70),
                dataLine.getDouble(ParameterTableColumn.POSTERIOR_80),
                dataLine.getDouble(ParameterTableColumn.POSTERIOR_90)));
        final PosteriorSummary posteriorSummary = new PosteriorSummary(center, lower, upper);
        posteriorSummary.setDeciles(deciles);
        return new AbstractMap.SimpleEntry<>(parameter, posteriorSummary);
    }
}
