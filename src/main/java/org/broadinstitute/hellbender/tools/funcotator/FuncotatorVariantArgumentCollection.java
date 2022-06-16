package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;

import java.io.File;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.stream.Collectors;

/**
 * Arguments to be be used by the {@link Funcotator} {@link org.broadinstitute.hellbender.engine.GATKTool},
 * which are specific to {@link Funcotator}.  Use this collection for small mutations (SNP, Indel)
 * Created by jonn on 9/12/18.
 */
public class FuncotatorVariantArgumentCollection extends BaseFuncotatorArgumentCollection {
    private static final long serialVersionUID = 1L;

    //-----------------------------------------------------
    // Required args:
    // (See superclass for more)


    //-----------------------------------------------------
    // Optional args:
    // (See superclass for more)

    @Argument(
            fullName = FuncotatorArgumentDefinitions.REMOVE_FILTERED_VARIANTS_LONG_NAME,
            optional = true,
            doc = "Ignore/drop variants that have been filtered in the input.  These variants will not appear in the output file."
    )
    public boolean removeFilteredVariants = false;

    @Argument(
            fullName = FuncotatorArgumentDefinitions.FIVE_PRIME_FLANK_SIZE_NAME,
            optional = true,
            doc = "Variants within this many bases of the 5' end of a transcript (and not overlapping any part of the transcript itself) will be annotated as being in the 5' flanking region of that transcript"
    )
    public int fivePrimeFlankSize = FuncotatorArgumentDefinitions.FIVE_PRIME_FLANK_SIZE_DEFAULT_VALUE;

    @Argument(
            fullName = FuncotatorArgumentDefinitions.THREE_PRIME_FLANK_SIZE_NAME,
            optional = true,
            doc = "Variants within this many bases of the 3' end of a transcript (and not overlapping any part of the transcript itself) will be annotated as being in the 3' flanking region of that transcript"
    )
    public int threePrimeFlankSize = FuncotatorArgumentDefinitions.THREE_PRIME_FLANK_SIZE_DEFAULT_VALUE;

    @Argument(
            fullName = FuncotatorArgumentDefinitions.REANNOTATE_VCF_LONG_NAME,
            optional = true,
            doc = "When input VCF has already been annotated, still annotate again."
    )
    public boolean reannotateVCF = false;
}
