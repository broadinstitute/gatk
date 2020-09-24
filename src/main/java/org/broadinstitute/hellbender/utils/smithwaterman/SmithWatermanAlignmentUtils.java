package org.broadinstitute.hellbender.utils.smithwaterman;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.AbstractReadThreadingGraph;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

/**
 * This class collects the various {@link SWParameters} that are used for various alignment procedures, along with
 * methods for reading them in from a TSV file.
 * It is likely that these parameters have not been rigorously optimized.
 * Documentation is also somewhat lacking, but we preserve any original comments that may have accompanied each set.
 * See also some relevant issues and comments:
 *  <p>
 *      <a href="http://github.com/broadinstitute/gatk/issues/2498">http://github.com/broadinstitute/gatk/issues/2498</a>
 *      <a href="http://github.com/broadinstitute/gatk/issues/5564">http://github.com/broadinstitute/gatk/issues/5564</a>
 *      <a href="http://github.com/broadinstitute/gatk/pull/4858#discussion_r194048530">http://github.com/broadinstitute/gatk/pull/4858#discussion_r194048530</a>
 *  </p>
 */
public final class SmithWatermanAlignmentUtils {
    /**
     * {@code ORIGINAL_DEFAULT} is only used in test code. It is worth noting that these tests are somewhat insensitive
     * to the particular values used; (e.g., the majority pass if {@link SmithWatermanAlignmentUtils#STANDARD_NGS}
     * is substituted, and all pass if {@link SmithWatermanAlignmentUtils#NEW_SW_PARAMETERS} is substituted).
     *
     * Original comments:
     *      match=1, mismatch = -1/3, gap=-(1+k/3)
     */
    public static final SWParameters ORIGINAL_DEFAULT = new SWParameters(3, -1, -4, -3);

    /**
     * {@code STANDARD_NGS} is the default for {@link AbstractReadThreadingGraph} methods for the recovery of dangling heads/tails.
     *
     * Original comments:
     *      none
     */
    public static final SWParameters STANDARD_NGS = new SWParameters(25, -50, -110, -6);

    /**
     * {@code NEW_SW_PARAMETERS} is the default for {@link CigarUtils#calculateCigar} for haplotype-to-reference alignment.
     * It was added in <a href="https://github.com/broadinstitute/gatk/pull/586">https://github.com/broadinstitute/gatk/pull/586</a>
     * at the same time as the {@link CigarUtils#calculateCigar} method. The original comments indicate that these values
     * were chosen via an optimization procedure, but no record or documentation of this procedure exists.
     * As indicated below in the original comments of {@link SmithWatermanAlignmentUtils#ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS},
     * the values chosen here heavily favor indels over substitutions.
     *
     * Original comments:
     *      used in the bubble state machine to apply Smith-Waterman to the bubble sequence
     *      these values were chosen via optimization against the NA12878 knowledge base
     */
    public static final SWParameters NEW_SW_PARAMETERS = new SWParameters(200, -150, -260, -11);

    /**
     * {@code ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS} is the default for read-to-haplotype alignment and was added in
     * <a href="https://github.com/broadinstitute/gatk/pull/4858">https://github.com/broadinstitute/gatk/pull/4858</a>
     * (superseding the use of {@link SmithWatermanAlignmentUtils#NEW_SW_PARAMETERS} in {@link AlignmentUtils#createReadAlignedToRef} in the
     * read-to-haplotype alignment step of that method).
     *
     * Original comments:
     *      In Mutect2 and HaplotypeCaller reads are realigned to their *best* haplotypes, which is very different from a generic alignment.
     *      The {@code NEW_SW_PARAMETERS} penalize a substitution error more than an indel up to a length of 9 bases!
     *      Suppose, for example, that a read has a single substitution error, say C -> T, on its last base.  Those parameters
     *      would prefer to extend a deletion until the next T on the reference is found in order to avoid the substitution, which is absurd.
     *      Since these parameters are for aligning a read to the biological sequence we believe it comes from, the parameters
     *      we choose should correspond to sequencer error.  They *do not* have anything to do with the prevalence of true variation!
     */
    public static final SWParameters ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS = new SWParameters(10, -15, -30, -5);

    @VisibleForTesting
    static final class SmithWatermanParametersTableReader extends TableReader<SWParameters> {
        @VisibleForTesting
        enum SmithWatermanParametersTableColumn {
            MATCH_VALUE("MATCH_VALUE"),
            MISMATCH_PENALTY("MISMATCH_PENALTY"),
            GAP_OPEN_PENALTY("GAP_OPEN_PENALTY"),
            GAP_EXTEND_PENALTY("GAP_EXTEND_PENALTY");

            private final String columnName;

            SmithWatermanParametersTableColumn(final String columnName) {
                this.columnName = Utils.nonNull(columnName);
            }

            @Override
            public String toString() {
                return columnName;
            }

            public static final TableColumnCollection COLUMNS = new TableColumnCollection(MATCH_VALUE, MISMATCH_PENALTY, GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY);
        }

        private SmithWatermanParametersTableReader(final Path path) throws IOException {
            super(path);
        }

        @Override
        protected SWParameters createRecord(final DataLine dataLine) {
            final int matchValue = dataLine.getInt(SmithWatermanParametersTableColumn.MATCH_VALUE);
            final int mismatchPenalty = dataLine.getInt(SmithWatermanParametersTableColumn.MISMATCH_PENALTY);
            final int gapOpenPenalty = dataLine.getInt(SmithWatermanParametersTableColumn.GAP_OPEN_PENALTY);
            final int gapExtendPenalty = dataLine.getInt(SmithWatermanParametersTableColumn.GAP_EXTEND_PENALTY);
            return new SWParameters(matchValue, mismatchPenalty, gapOpenPenalty, gapExtendPenalty);
        }
    }

    /**
     * Reads in {@link SWParameters} from a TSV file. The file must contain a single column-header line corresponding
     * to the columns defined in {@link SmithWatermanParametersTableReader.SmithWatermanParametersTableColumn}
     * and a single row of integer parameter values.
     */
    public static SWParameters readSmithWatermanParametersFromTSV(final GATKPath path) {
        Utils.nonNull(path);
        try (final SmithWatermanParametersTableReader tableReader = new SmithWatermanParametersTableReader(path.toPath())) {
            final List<SWParameters> swParametersList = tableReader.toList();
            if (swParametersList.size() != 1) {
                throw new UserException.BadInput(String.format(
                        "TSV file should contain a single column-header row and a single row of Smith-Waterman parameter values: %s", path.toPath()));
            }
            return swParametersList.get(0);
        } catch (final IOException e) {
            throw new UserException.BadInput(String.format("Could not read Smith-Waterman parameters from TSV file: %s", path.toPath()));
        }
    }
}
