package org.broadinstitute.hellbender.transformers;

import org.apache.logging.log4j.util.Supplier;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * A read transformer to convert IUPAC bases (i.e. non-ATCGs) to Ns
 * Some references (like human hg38) contain IUPAC bases that can be propagated into the reads when decoding cram
 * This transformation is done in-place without copying
 */
public class IUPACReadTransformer implements ReadTransformer {
    private static final long serialVersionUID = 1L;
    private boolean strictMode;
    private OneShotLogger logger = new OneShotLogger(this.getClass());

    public IUPACReadTransformer() {
        this.strictMode = false;
    }

    public IUPACReadTransformer(final boolean strictMode) {
        this.strictMode = strictMode;
    }

    @Override
    public GATKRead apply(GATKRead read) {
        final byte[] maybeTransformed = BaseUtils.convertIUPACtoN(read.getBases(), strictMode, false);
        if (!Arrays.equals(read.getBases(), maybeTransformed)) {
            logger.warn(() -> "At least one read contains IUPAC bases that have been transformed.  Read " + read.getName() + " contains: "
            + IntStream.range(0, read.getBases().length).map(idx -> read.getBase(idx))
                    .filter(i -> !BaseUtils.isNucleotide((byte)i) && !BaseUtils.isNBase((byte)i))
                    .mapToObj(i -> (char)i).collect(Collectors.toList()));
            read.setBases(maybeTransformed);
        }
        return read;
    }
}
