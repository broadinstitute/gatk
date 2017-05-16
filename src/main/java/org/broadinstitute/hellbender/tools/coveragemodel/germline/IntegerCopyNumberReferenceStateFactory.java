package org.broadinstitute.hellbender.tools.coveragemodel.germline;

import org.broadinstitute.barclay.utils.Utils;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.IntegerCopyNumberState;
import org.broadinstitute.hellbender.tools.exome.sexgenotyper.GermlinePloidyAnnotatedTargetCollection;
import org.broadinstitute.hellbender.tools.exome.sexgenotyper.SexGenotypeData;

import javax.annotation.Nonnull;
import java.io.Serializable;
import java.util.function.BiFunction;

/**
 * This class implements a {@link BiFunction} that takes a sex genotype and a target and outputs the reference
 * integer copy number state, as an instance of {@link IntegerCopyNumberState}, for that sex genotype and target.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class IntegerCopyNumberReferenceStateFactory implements BiFunction<SexGenotypeData, Target,
        IntegerCopyNumberState>, Serializable {

    private static final long serialVersionUID = 6322699892999789289L;

    private final GermlinePloidyAnnotatedTargetCollection germlinePloidyAnnotatedTargetCollection;

    public IntegerCopyNumberReferenceStateFactory(@Nonnull final GermlinePloidyAnnotatedTargetCollection
                                                          germlinePloidyAnnotatedTargetCollection) {
        this.germlinePloidyAnnotatedTargetCollection = Utils.nonNull(germlinePloidyAnnotatedTargetCollection);
    }

    @Override
    public IntegerCopyNumberState apply(@Nonnull final SexGenotypeData sexGenotypeData, @Nonnull final Target target) {
        return new IntegerCopyNumberState(germlinePloidyAnnotatedTargetCollection
                .getTargetGermlinePloidyByGenotype(target, sexGenotypeData.getSexGenotype()));
    }
}
