package org.broadinstitute.hellbender.tools.spark.longread;


import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.LongReadAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;

import java.util.Arrays;
import java.util.List;

/**
 * Amend CIGARs of long reads showing the following strange behaviors:
 *
 * <p>
 *     Currently we define "strange" as:
 *     <ul>
 *         <li>CIGAR having clipping neighboring insertion blocks (alignments will have a tag "YS:i:1"),
 *             CIGARs of this type will be amended in a way such that the bases marked as insertion will be marked
 *             as soft-clipped (hence be merged with the already soft-clipped bases, if that applies)
 *         </li>
 *         <li>CIGAR having an insertion/deletion operation neighboring a deletion/insertion operation
 *             (alignments will have a tag "YG:Z:section_of_strange_cigar", and a tag "YO:Z:original_cigar")
 *             CIGARs of this type will be amended in a way such that original a(I|D)b(D|I) will become
 *             abs(a-b)(I|D)min(a,b)M.
 *             Very importantly,
 *             <ul>
 *                 <li>the MD tag will be removed as the original value will be invalidated
 *                     after such operation, and re-computing them will be impractical.
 *                 </li>
 *                 <li>the AS tag will NOT be modified; while technically this is wrong, we know some tools
 *                     downstream need this tag,
 *                     and we HOPE the original value will not be so far away from the "good" value.
 *                 </li>
 *             </ul>
 *
 *
 *         </li>
 *     </ul>
 * </p>
 */
@DocumentedFeature
@CommandLineProgramProperties(summary = "Amend CIGARs of alignment records that seem strange.",
        oneLineSummary = "Amend CIGARs of alignment records that seem strange",
        programGroup = LongReadAnalysisProgramGroup.class)
@BetaFeature
public class SimpleAmendLongReadsCigarSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    public static final String NEIGHBOR_CLIPGAP_OP = "YS";
    public static final String NEIGHBOR_INSDEL_OP = "YG";
    public static final String NEIGHBOR_INSDEL_ORIGINAL_CIGAR = "YO";


    @Override
    public boolean requiresReads()
    {
        return true;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Arrays.asList(new ReadFilterLibrary.AllowAllReadsReadFilter());
    }

    @Override
    protected void runTool( final JavaSparkContext ctx ) {

    }
}
