package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.collections4.list.SetUniqueList;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.hdf5.PoN;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Utility class use to map matrices rows between two different target lists.
 * <p>
 *     The "case" target list is provided by a {@link ReadCountCollection}, whereas the "pon" target list
 *     is provided by a {@link PoN}.
 * </p>
 * <p>
 *     This class contains methods to re-organize matrix rows, based on one list, to match the order
 *     imposed by the other list.
 * </p>
 * <p>
 *     It also provides methods to create a PoN re-organized count matrix directly from a {@link ReadCountCollection}.
 *     object and the reverse operation, that is to create a new {@link ReadCountCollection} from a count matrix
 *     whose rows have been reorganized to satisfy the {@link PoN} target order.
 * </p>
 */
public final class Case2PoNTargetMapper {

    private final List<String> ponTargetNames;
    private final Map<String, Integer> caseTargetIndexes;
    private final List<Target> outputTargets;

    /**
     * Creates a new mapper.
     * @param caseTargets the case sample targets as they appear in the input read counts.
     * @param ponTargetNames the PoN target names in the order they are present in the PoN.
     */
    public Case2PoNTargetMapper(final List<Target> caseTargets, final List<String> ponTargetNames)  {

        this.ponTargetNames = ponTargetNames;
        this.caseTargetIndexes = IntStream.range(0, caseTargets.size()).boxed()
                .collect(Collectors.toMap(i -> caseTargets.get(i).getName(), Function.identity()));
        if (this.ponTargetNames.stream().anyMatch(n -> ! caseTargetIndexes.containsKey(n))) {
            throw new UserException.BadInput(String.format("the read-count input file is missing some target in the PoN: e.g. '%s'",
                    this.ponTargetNames.stream().filter(n -> ! caseTargetIndexes.containsKey(n)).limit(10).collect(Collectors.joining(", "))));
        }
        this.outputTargets =  ponTargetNames.stream()
                .map(caseTargetIndexes::get).sorted().map(caseTargets::get).collect(Collectors.toList());
    }

    /**
     * Re-arrange case read counts counts based on PoN target order.
     * @param caseCounts the input counts to rearrange.
     * @return never {@code null}, a new matrix with row sorted according to the PoN target order.
     */
    public RealMatrix fromCaseToPoNCounts(final RealMatrix caseCounts) {
        final RealMatrix result = new Array2DRowRealMatrix(ponTargetNames.size(), caseCounts.getColumnDimension());
        for (int i = 0; i < ponTargetNames.size(); i++) {
            final String targetName = ponTargetNames.get(i);
            final Integer inputIndex = caseTargetIndexes.get(targetName);
            if (inputIndex == null) {
                throw new UserException.BadInput(String.format("missing PoN target name %s in input read counts", targetName));
            }
            result.setRow(i, caseCounts.getRow(inputIndex));
        }
        return result;
    }

    /**
     * Re-arrange the input rows from the PoN to the case data target order.
     * @param ponCounts count matrix with row organized using the PoN target order.
     * @return never {@code null} a new matrix with the row order changed according to the case read count target order.
     */
    public RealMatrix fromPoNToCaseCounts(final RealMatrix ponCounts) {
        final Map<String,Integer> ponTargetsIndexes = IntStream.range(0, ponTargetNames.size()).boxed()
                .collect(Collectors.toMap(ponTargetNames::get, Function.identity()));
        final RealMatrix result = new Array2DRowRealMatrix(ponCounts.getRowDimension(),ponCounts.getColumnDimension());
        for (int i = 0; i < outputTargets.size(); i++) {
            final Target target = outputTargets.get(i);
            result.setRow(i, ponCounts.getRow(ponTargetsIndexes.get(target.getName())));
        }
        return result;
    }

    /**
     * Given a read counts using the PoN target order to a read-count collection
     * based on the case sample target order.
     *
     * @param counts PoN target sorted count matrix.
     * @param countColumnNames count column names.
     * @return never {@code null} a read-count collection with the row order according to the case read count target order.
     */
    public ReadCountCollection fromPoNtoCaseCountCollection(final RealMatrix counts, final List<String> countColumnNames) {
        return new ReadCountCollection(SetUniqueList.setUniqueList(new ArrayList<>(outputTargets)),SetUniqueList.setUniqueList(new ArrayList<>(countColumnNames)),
                fromPoNToCaseCounts(counts));
    }
}