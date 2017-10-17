package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.tools.walkers.readorientation.AltSiteRecord.AltSiteRecordTableWriter;
import org.broadinstitute.hellbender.tools.walkers.readorientation.RefSiteHistogram.RefSiteHistogramWriter;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by Takuto Sato on 7/26/17.
 */

@CommandLineProgramProperties(
        summary = "Traverse and collect data from a bam for Mutect2 Read Orientation Filter",
        oneLineSummary = "Data collection for Mutect2 Read Orientation Filter",
        programGroup = VariantProgramGroup.class
)

public class CollectDataForReadOrientationFilter extends LocusWalker {
    public static final String ALT_DATA_TABLE_SHORT_NAME = "alt_table";
    public static final String ALT_DATA_TABLE_LONG_NAME = "alt_data_table";

    public static final String REF_HISTOGRAM_TABLE_SHORT_NAME = "ref_table";
    public static final String REF_HISTOGRAM_TABLE_LONG_NAME = "ref_histogram_table";


    @Argument(fullName = "min-mq", shortName = "mq", doc = "exclude reads below this quality from pileup", optional = true)
    static int MINIMUM_MEDIAN_MQ = 20;

    @Argument(fullName = "min-bq", shortName = "bq", doc = "exclude bases below this quality from pileup", optional = true)
    static int MINIMUM_BASE_QUALITY = 10;

    @Argument(fullName = ALT_DATA_TABLE_LONG_NAME,
            shortName = ALT_DATA_TABLE_SHORT_NAME,
            doc = "a tab-separated output table of pileup data over alt sites")
    static File altDataTable = null;

    @Argument(fullName = REF_HISTOGRAM_TABLE_LONG_NAME,
            shortName = REF_HISTOGRAM_TABLE_SHORT_NAME,
            doc = "a tab-separated output histogram of depths over ref sites")
    static File refHistogramTable = null;

    // for computational efficiency, for each reference context, we build a depth histogram over ref sites
    private static RefSiteHistogram[] refSiteHistograms = new RefSiteHistogram[64];

    AltSiteRecordTableWriter altTableWriter;

    @Override
    public boolean requiresReference(){
        return true;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return ReadUtils.makeStandardReadFilters();
    }

    @Override
    public void onTraversalStart() {
        final List<String> all3Mers = SequenceUtil.generateAllKmers(3).stream().map(String::new).collect(Collectors.toList());
        for (final String refContext : all3Mers){
            refSiteHistograms[contextToIndex(refContext)] = new RefSiteHistogram(refContext);
        }

        // intentionally not use try-with-resources so that the writer stays open outside of the try block
        try {
            altTableWriter = new AltSiteRecordTableWriter(altDataTable);
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception creating a writer for %s", altDataTable), e);
        }

    }

    @Override
    public void apply(final AlignmentContext alignmentContext, final ReferenceContext referenceContext, final FeatureContext featureContext){
        // referenceContext always comes with a window of a single base, so
        // manually expand the window and get the 3-mer for now.
        // TODO: implement getBasesInInterval() in referenceContext. Maybe simplify to getKmer(int k)?
        // TODO: this is still relevant (10/2). I shouldn't mess with the internal state of the ref context object
        referenceContext.setWindow(1, 1);
        final String refContext = new String(referenceContext.getBases());
        assert refContext.length() == 3 : "kmer must have length 3";
        if (refContext.contains("N")) {
            return;
        }

        if (refContext == null){
            logger.info(String.format("null reference found at interval %s, k-mer = %s",
                    referenceContext.getInterval().toString(), refContext));
            return;
        }

        final ReadPileup pileup = alignmentContext.filter(pe -> pe.getQual() > MINIMUM_BASE_QUALITY).getBasePileup();

        // This case should not happen, as AlignmentContext should come filtered, it does happen once in a while
        if (pileup.size() == 0){
            return;
        }

        /*** A series of heuristics to detect a bad pileup ***/
        // skip indels
        if (pileup.getNumberOfElements(pe -> pe.isDeletion() || pe.isAfterInsertion() || pe.isBeforeDeletionStart()) > 0){
            return;
        }

        // skip low MQ loci
        List<Integer> mappingQualities = new ArrayList<>(pileup.size());

        // There is no shortcut or a standard API for converting an int[] to List<Integer> (Arrays.toList gives us List<int[]>)
        // so we must do so manually with a for loop
        // Median in Apache commons would be an alternative approach, but it takes in a double[], not int[],
        // so converting to double[] will be just as expensitve as iterating
        for (final int mq : pileup.getMappingQuals()) {
            mappingQualities.add(mq);
        }

        final double medianMQ = MathUtils.median(mappingQualities);

        if (medianMQ < MINIMUM_MEDIAN_MQ) {
            return;
        }

        final int[] baseCounts = pileup.getBaseCounts();
        final int depth = (int) MathUtils.sum(baseCounts);

        if (depth == 0) {
            if (pileup.size() != 0) {
                logger.info(String.format("Encountered depth = 0 at %s. Pileup size was %d", referenceContext.getInterval(), pileup.size()));
            }
            return;
        }

        final Nucleotide refBase = Nucleotide.valueOf(refContext.getBytes()[1]);

        // Make a copy of base counts and update the counts of ref to -1. Now the argmax of the array gives us
        // the alt base.
        final int[] baseCountsCopy = Arrays.copyOf(baseCounts, baseCounts.length);
        baseCountsCopy[refBase.ordinal()] = -1;
        final int altBaseIndex = MathUtils.argmax(baseCountsCopy);
        final boolean referenceSite = baseCounts[altBaseIndex] == 0;

        /*** End heuristics ***/

        // if the site is ref, we simply update the coverage histogram
        if (referenceSite){
            refSiteHistograms[contextToIndex(refContext)].increment(depth);
            return;
        }

        // if we got here, we have an alt site
        final Nucleotide altBase = Nucleotide.valueOf(BaseUtils.baseIndexToSimpleBase(altBaseIndex));

        final int[] altF1R2Counts = new int[4];
        for (Nucleotide base : Arrays.asList(Nucleotide.A, Nucleotide.C, Nucleotide.G, Nucleotide.G)){
            altF1R2Counts[base.ordinal()] = pileup.getNumberOfElements(pe ->
                    Nucleotide.valueOf(pe.getBase()) == base && ! ReadUtils.isF2R1(pe.getRead()));
        }

        try {
            altTableWriter.writeRecord(new AltSiteRecord(refContext, baseCounts, altF1R2Counts, depth, altBase));
        } catch (IOException e) {
            throw new UserException("Encountered an IO Exception writing to the alt data table", e);
        }

        return;
    }

    @Override
    public Object onTraversalSuccess() {
        RefSiteHistogram.writeRefSiteHistograms(Arrays.asList(refSiteHistograms), refHistogramTable);
        return "SUCCESS";
    }

    @Override
    public void closeTool() {
        if (altTableWriter != null) {
            try {
                altTableWriter.close();
            } catch (IOException e) {
                throw new UserException("Encountered an IO exception while closing the alt table writer", e);
            }
        }
    }


    /***
     * Maps a reference 3-mer to an array index using the quaternary (base-4) numeral system
     * Example: AGT is represented as 023 in quaternary, which in decimal is 0*4^2 + 2*4^1 + 3*4^0 = 0+8+3 = 11
     *
     * @param reference3mer
     * @return
     */
    protected int contextToIndex(String reference3mer){
        Utils.validateArg(reference3mer.matches("[ACGT]{3}"),
                "input must be a string of length 3 and comprise of A, C, G,and T");
        final int digit2 = Nucleotide.valueOf(reference3mer.substring(0, 1)).ordinal();
        final int digit1 = Nucleotide.valueOf(reference3mer.substring(1, 2)).ordinal();
        final int digit0 = Nucleotide.valueOf(reference3mer.substring(2, 3)).ordinal();

        return digit2 * 16 + digit1 * 4 + digit0;
    }
}
