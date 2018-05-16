package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordComparator;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import scala.Tuple2;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.Serializable;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.InversionBreakendPreFilter.OverlappingPair;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.SYMB_ALT_ALLELE_INV;

/**
 * Depending on how the two intervals overlap each other (see {@link OverlappingScenario}),
 * the signatures tell us what possible types of events are involved.
 * However, there are possibly some degeneracies in some of the overlap scenarios
 * that needs to be resolved.
 */
final class LinkedInversionBreakpointsInference implements Serializable {
    private static final long serialVersionUID = 1L;

    public enum OverlappingScenario {
        // what it means for possible variation(s) is described in toString()
        THREE_CONTAINS_FIVE,      // interval from INV33 contains that from INV55
        FIVE_CONTAINS_THREE,      // interval from INV55 contains that from INV33
        THREE_INTERLEAVES_FIVE,   // interval from INV33 upstream of that from INV55, but doesn't contain it
        FIVE_INTERLEAVES_THREE;   // interval from INV55 upstream of that from INV33, but doesn't contain it

        @Override
        public String toString() { // a "primed-annotated" block is inverted
            switch (this) {
                case THREE_CONTAINS_FIVE:
                    return " ABC -> A(B|B')A' ";
                case FIVE_CONTAINS_THREE:
                    return " ABC -> C'(B|B')C ";
                case THREE_INTERLEAVES_FIVE:
                    return " ABC -> AC'B'A'C  ";
                case FIVE_INTERLEAVES_THREE:
                    return " ABC ->    B'     ";
                default: throw new IllegalStateException();
            }
        }

        String getDescription() {
            return toString();
        }
    }

    //==================================================================================================================

    /**
     * See {@link OverlappingScenario}
     */
    static List<VariantContext> makeInterpretation(final List<OverlappingPair> overlappingPairs,
                                                   final String pathToFastqsDir,
                                                   final ReferenceMultiSource reference,
                                                   final Logger toolLogger) {

        return Collections.emptyList();
    }
}
