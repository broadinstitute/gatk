package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleMetadata;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.Decile;
import org.broadinstitute.hellbender.utils.mcmc.DecileCollection;
import org.broadinstitute.hellbender.utils.mcmc.ParameterEnum;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class ParameterDecileCollection<T extends Enum<T> & ParameterEnum> extends AbstractSampleRecordCollection<Map.Entry<T, DecileCollection>> {
    //note to developers: repeat the column headers in Javadoc so that they are viewable when linked
    /**
     * PARAMETER_NAME, POSTERIOR_10, POSTERIOR_20, POSTERIOR_30, POSTERIOR_40, POSTERIOR_50, POSTERIOR_60, POSTERIOR_70, POSTERIOR_80, POSTERIOR_90
     */
    enum ParameterTableColumn {
        PARAMETER_NAME,
        POSTERIOR_10,
        POSTERIOR_20,
        POSTERIOR_30,
        POSTERIOR_40,
        POSTERIOR_50,
        POSTERIOR_60,
        POSTERIOR_70,
        POSTERIOR_80,
        POSTERIOR_90;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static DecileCollection parseDecilesFromDataLine(final DataLine dataLine) {
        return new DecileCollection(Arrays.asList(
                dataLine.getDouble(ParameterTableColumn.POSTERIOR_10),
                dataLine.getDouble(ParameterTableColumn.POSTERIOR_20),
                dataLine.getDouble(ParameterTableColumn.POSTERIOR_30),
                dataLine.getDouble(ParameterTableColumn.POSTERIOR_40),
                dataLine.getDouble(ParameterTableColumn.POSTERIOR_50),
                dataLine.getDouble(ParameterTableColumn.POSTERIOR_60),
                dataLine.getDouble(ParameterTableColumn.POSTERIOR_70),
                dataLine.getDouble(ParameterTableColumn.POSTERIOR_80),
                dataLine.getDouble(ParameterTableColumn.POSTERIOR_90)));
    }

    private static void appendDecilesToDataLine(final DataLine dataLine,
                                                final DecileCollection deciles) {
        dataLine.append(formatDouble(deciles.get(Decile.DECILE_10)))
                .append(formatDouble(deciles.get(Decile.DECILE_20)))
                .append(formatDouble(deciles.get(Decile.DECILE_30)))
                .append(formatDouble(deciles.get(Decile.DECILE_40)))
                .append(formatDouble(deciles.get(Decile.DECILE_50)))
                .append(formatDouble(deciles.get(Decile.DECILE_60)))
                .append(formatDouble(deciles.get(Decile.DECILE_70)))
                .append(formatDouble(deciles.get(Decile.DECILE_80)))
                .append(formatDouble(deciles.get(Decile.DECILE_90)));
    }

    private final Map<T, DecileCollection> parameterToDecileCollectionMap;

    public ParameterDecileCollection(final SampleMetadata metadata,
                                     final Map<T, DecileCollection> parameterToDecileCollectionMap,
                                     final Class<T> parameterClass) {
        super(
                Utils.nonNull(metadata),
                new ArrayList<>(parameterToDecileCollectionMap.entrySet()),
                ParameterTableColumn.COLUMNS,
                dataLine -> {
                    final String parameterName = dataLine.get(ParameterTableColumn.PARAMETER_NAME);
                    final T parameter = Enum.valueOf(Utils.nonNull(parameterClass), parameterName);
                    final DecileCollection deciles = parseDecilesFromDataLine(dataLine);
                    return new AbstractMap.SimpleEntry<>(parameter, deciles);},
                (record, dataLine) -> {
                    final T parameter = record.getKey();
                    final DecileCollection deciles = record.getValue();
                    appendDecilesToDataLine(dataLine.append(parameter.toString()), deciles);
                }
        );
        this.parameterToDecileCollectionMap = parameterToDecileCollectionMap;
    }

    public ParameterDecileCollection(final File file,
                                     final Class<T> parameterClass) {
        super(
                Utils.nonNull(file),
                ParameterTableColumn.COLUMNS,
                dataLine -> {
                    final String parameterName = dataLine.get(ParameterTableColumn.PARAMETER_NAME);
                    final T parameter = Enum.valueOf(Utils.nonNull(parameterClass), parameterName);
                    final DecileCollection deciles = parseDecilesFromDataLine(dataLine);
                    return new AbstractMap.SimpleEntry<>(parameter, deciles);},
                (record, dataLine) -> {
                    final T parameter = record.getKey();
                    final DecileCollection deciles = record.getValue();
                    dataLine.append(parameter.toString());
                    appendDecilesToDataLine(dataLine, deciles);
                }
        );
        parameterToDecileCollectionMap = getRecords().stream().collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
    }

    public DecileCollection getDeciles(final T parameter) {
        return parameterToDecileCollectionMap.get(parameter);
    }
}