package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPathSpecifier;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;

import java.io.IOException;
import java.util.*;

/**
 * Example/toy program that shows how to implement the VariantWalker interface. Prints supplied variants
 * along with overlapping reads/reference bases/variants (if present).
 */
@CommandLineProgramProperties(
        summary = "Example tool that prints variants supplied to the specified output file (stdout if none provided), along with overlapping reads/reference bases/variants (if provided)",
        oneLineSummary = "Example tool that prints variants with optional contextual data",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public final class BlahVariantWalker extends VariantWalker {
    static final Logger logger = LogManager.getLogger(BlahVariantWalker.class);

    private final int GQ_CUTOFF = 10000;
    private final char SEPARATOR = '\t';
    private SimpleXSVWriter vetWriter = null;
    private SimpleXSVWriter petWriter = null;

    private GenomeLocSortedSet intervalArgumentGenomeLocSortedSet;
    private GenomeLocSortedSet coverageLocSortedSet;
    private SimpleInterval previousInterval;
    private String sampleName;


    @Argument(fullName = "vet-table-out-path",
            shortName = "VO",
            doc="Path to where the variants expanded table should be written")
    public GATKPathSpecifier vetOutput = null;


    @Argument(fullName = "pet-table-out-path",
            shortName = "PO",
            doc="Path to where the positions table should be written")
    public GATKPathSpecifier petOutput = null;

    public BlahPetCreation.GQStateEnum gqStateToIgnore = null;

    @Override
    public boolean requiresIntervals() {
        return true;
    }

    @Override
    public void onTraversalStart() {

        final SAMSequenceDictionary seqDictionary = getBestAvailableSequenceDictionary();

        final GenomeLocSortedSet genomeLocSortedSet = new GenomeLocSortedSet(new GenomeLocParser(seqDictionary));
        intervalArgumentGenomeLocSortedSet = GenomeLocSortedSet.createSetFromList(genomeLocSortedSet.getGenomeLocParser(), IntervalUtils.genomeLocsFromLocatables(genomeLocSortedSet.getGenomeLocParser(), intervalArgumentCollection.getIntervals(seqDictionary)));
        coverageLocSortedSet = new GenomeLocSortedSet(new GenomeLocParser(seqDictionary));

        try {
            List<String> vetHeader = BlahVetCreation.getHeaders();
            vetWriter = new SimpleXSVWriter(vetOutput.toPath(), SEPARATOR);
            vetWriter.setHeaderLine(vetHeader);

        } catch (final IOException e) {
            throw new IllegalArgumentException("Current variant is missing required fields", e);
        }
        try {
            List<String> petHeader = BlahPetCreation.getHeaders();
            petWriter = new SimpleXSVWriter(petOutput.toPath(), SEPARATOR);
            petWriter.setHeaderLine(petHeader);

        } catch (final IOException e) {
            throw new IllegalArgumentException("Current variant is missing required fields", e);
        }
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        if (sampleName == null) {
            sampleName = variant.getGenotype(0).getSampleName();
        }

        // get the intervals this variant covers
        final GenomeLoc variantGenomeLoc = intervalArgumentGenomeLocSortedSet.getGenomeLocParser().createGenomeLoc(variant.getContig(), variant.getStart(), variant.getEnd());
        final List<GenomeLoc> intervalsToWrite = intervalArgumentGenomeLocSortedSet.getOverlapping(variantGenomeLoc);

        if (intervalsToWrite.size() == 0){
            throw new IllegalStateException("There are no intervals being covered by this variant");
        }

        // take the first interval(assuming this is returned in order) and make sure if its a variant, that it starts at/after the interval start
        // we are going to ignore any deletions that start before an interval.
        if (!variant.isReferenceBlock() && intervalsToWrite.get(0).getStart() > variant.getStart()){
            return;
        }

        // if the only alt allele for a variant is `*`, we ignore it
        if (!variant.isReferenceBlock() &&  variant.getAlternateAlleles().size() == 2 && variant.hasAlternateAllele(Allele.SPAN_DEL)){
            return;
        }

        // create VET output
        if (!variant.isReferenceBlock()) {
            final List<String> TSVLineToCreateVet = BlahVetCreation.createVariantRow(variant);

            // write the variant to the XSV
            SimpleXSVWriter.LineBuilder vetLine = vetWriter.getNewLineBuilder();
            vetLine.setRow(TSVLineToCreateVet);
            vetLine.write();
        }
        boolean firstInterval = true;
        for (GenomeLoc genomeLoc : intervalsToWrite) {

            int start = Math.max(genomeLoc.getStart(), variant.getStart());
            int end = Math.min(genomeLoc.getEnd(), variant.getEnd());

            // create PET output if the reference block's GQ is not the one to throw away or its a variant
            if (!variant.isReferenceBlock() || !BlahPetCreation.getGQStateEnum(variant.getGenotype(0).getGQ()).equals(gqStateToIgnore)) {

                // add interval to "covered" intervals
                // GenomeLocSortedSet will automatically merge intervals that are overlapping when setting `mergeIfIntervalOverlaps`
                // to true.  In a GVCF most blocks are adjacent to each other so they wouldn't normally get merged.  We check
                // if the current record is adjacent to the previous record and "overlap" them if they are so our set is as
                // small as possible while still containing the same bases.
                final SimpleInterval variantInterval = new SimpleInterval(variant.getContig(), start, end);
                final int intervalStart = (previousInterval != null && previousInterval.overlapsWithMargin(variantInterval, 1)) ?
                        previousInterval.getStart() : variantInterval.getStart();
                final int intervalEnd = (previousInterval != null && previousInterval.overlapsWithMargin(variantInterval, 1)) ?
                        Math.max(previousInterval.getEnd(), variantInterval.getEnd()) : variantInterval.getEnd();
                final GenomeLoc possiblyMergedGenomeLoc = coverageLocSortedSet.getGenomeLocParser().createGenomeLoc(variantInterval.getContig(), intervalStart, intervalEnd);
                coverageLocSortedSet.add(possiblyMergedGenomeLoc, true);
                previousInterval = new SimpleInterval(possiblyMergedGenomeLoc);

                List<List<String>> TSVLinesToCreatePet;
                // handle deletions that span across multiple intervals
                if (!firstInterval && !variant.isReferenceBlock()) {
                    TSVLinesToCreatePet = BlahPetCreation.createSpanDelRows(start, variant, end);
                } else {
                    TSVLinesToCreatePet = BlahPetCreation.createPositionRows(start, variant, end);
                }

                // write the position to the XSV
                for (List<String> TSVLineToCreatePet : TSVLinesToCreatePet) {
                    petWriter.getNewLineBuilder().setRow(TSVLineToCreatePet).write();
                }
            }
            firstInterval = false;
        }
    }

    @Override
    public Object onTraversalSuccess() {
        final GenomeLocSortedSet uncoveredIntervals = intervalArgumentGenomeLocSortedSet.subtractRegions(coverageLocSortedSet);
        logger.info(uncoveredIntervals.coveredSize());
        logger.info((1.0 * uncoveredIntervals.coveredSize()) / intervalArgumentGenomeLocSortedSet.coveredSize() );
        for (GenomeLoc genomeLoc : uncoveredIntervals) {
            // write the position to the XSV
            for (List<String> TSVLineToCreatePet : BlahPetCreation.createMissingTSV(genomeLoc.getStart(), genomeLoc.getEnd(),sampleName)) {
                petWriter.getNewLineBuilder().setRow(TSVLineToCreatePet).write();
            }
        }
        return 0;
    }

    @Override
    public void closeTool() {
        if (vetWriter != null ) {
            try {
                vetWriter.close();
            } catch (final Exception e) {
                throw new IllegalArgumentException("Couldn't close VET writer", e);
            }
        }
        if (petWriter != null) {
            try {
                petWriter.close();
            } catch (final Exception e) {
                throw new IllegalArgumentException("Couldn't close PET writer", e);
            }
        }
    }
}
