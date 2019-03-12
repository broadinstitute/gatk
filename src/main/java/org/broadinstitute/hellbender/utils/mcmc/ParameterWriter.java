package org.broadinstitute.hellbender.utils.mcmc;

import java.nio.file.Path;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.Map;

/**
 * Table writer for ouputting a {@link PosteriorSummary} for each global parameter in a {@link ParameterEnum}.
 *
 * @author Samuel Lee &lt;valentin@broadinstitute.org&gt;
 */
public final class ParameterWriter<T extends Enum<T> & ParameterEnum> extends TableWriter<Map.Entry<T, PosteriorSummary>> {
    
    private final String doubleFormat;

    public ParameterWriter(final Path path, final String doubleFormat) throws IOException {
        super(path, ParameterTableColumn.COLUMNS);
        this.doubleFormat = doubleFormat;
    }

    @Override
    protected void composeLine(final Map.Entry<T, PosteriorSummary> record, final DataLine dataLine) {
        final T parameter = record.getKey();
        final PosteriorSummary posteriorSummary = record.getValue();
        dataLine.append(parameter.toString())
                .append(formatDouble(posteriorSummary.getCenter()))
                .append(formatDouble(posteriorSummary.getLower()))
                .append(formatDouble(posteriorSummary.getUpper()))
                .append(formatDouble(posteriorSummary.getDeciles().get(Decile.DECILE_10)))
                .append(formatDouble(posteriorSummary.getDeciles().get(Decile.DECILE_20)))
                .append(formatDouble(posteriorSummary.getDeciles().get(Decile.DECILE_30)))
                .append(formatDouble(posteriorSummary.getDeciles().get(Decile.DECILE_40)))
                .append(formatDouble(posteriorSummary.getDeciles().get(Decile.DECILE_50)))
                .append(formatDouble(posteriorSummary.getDeciles().get(Decile.DECILE_60)))
                .append(formatDouble(posteriorSummary.getDeciles().get(Decile.DECILE_70)))
                .append(formatDouble(posteriorSummary.getDeciles().get(Decile.DECILE_80)))
                .append(formatDouble(posteriorSummary.getDeciles().get(Decile.DECILE_90)));
    }

    private String formatDouble(final double d) {
        return String.format(doubleFormat, d);
    }
}
