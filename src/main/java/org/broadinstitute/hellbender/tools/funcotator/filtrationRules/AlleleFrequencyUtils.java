package org.broadinstitute.hellbender.tools.funcotator.filtrationRules;

import org.broadinstitute.hellbender.tools.funcotator.FilterFuncotations.AlleleFrequencyDataSource;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Allele frequency utilities that are dataset-agnostic
 */
public abstract class AlleleFrequencyUtils {

    /**
     * Build a {@link FuncotationFiltrationRule} matching Funcotations from variants with a
     * maximum MAF less than some threshold.
     *
     * @param maxMaf the MAF threshold to check in the rule. Must be in the range [0, 1]
     * @param afDataSource the allele frequency data source (ExAC or gnomAD) with which the original VCF was Funcotated.
     * @return a {@link FuncotationFiltrationRule} matching Funcotations with a MAF
     *         less than {@code maxMaf} across all sub-populations.
     */
    public static FuncotationFiltrationRule buildMaxMafRule(final double maxMaf, final AlleleFrequencyDataSource afDataSource) {
        ParamUtils.inRange(maxMaf, 0, 1, "MAF must be between 0 and 1");
        if (afDataSource.equals(AlleleFrequencyDataSource.exac)) {
            return (funcotations, variant) -> AlleleFrequencyExacUtils.getMaxMinorAlleleFreq(funcotations) <= maxMaf;
        }
        else {
            return (funcotations, variant) -> {
                final Map<String, List<String>> condensedFuncotations = funcotations.stream().collect(
                        Collectors.groupingBy(Map.Entry::getKey,
                                              Collectors.mapping(Map.Entry::getValue, Collectors.toList()))
                );
                return (!AlleleFrequencyGnomadUtils.allFrequenciesFiltered(condensedFuncotations)
                        && AlleleFrequencyGnomadUtils.getMaxMinorAlleleFreq(condensedFuncotations) <= maxMaf);
            };
        }
    }
}
