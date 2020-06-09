package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.util.ArrayList;
import java.util.List;

/**
 * Replace bases in reads with reference bases.
 *
 * Used to sanitize data from samples for the purposes of eliminating Personal Identifiable Information.
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        summary = "Replace bases in reads with reference bases.",
        oneLineSummary = "Replace bases in reads with reference bases.",
        programGroup = OtherProgramGroup.class
)
public final class ReadSanitizer extends ReadWalker {

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Output bam file.")
    public GATKPath OUTPUT;

    @Argument(
            fullName = "ref-base-quality",
            shortName = "ref-base-quality",
            doc = "Quality for bases that are set to the reference base.",
            minValue = 0,
            maxValue = 60,
            optional = true
    )
    public int refQual = 60;

    @Argument(
            fullName = "use-extended-cigar",
            shortName = "use-extended-cigar",
            doc = "If true, will produce `=` instead of `M` for matching bases.",
            optional = true
    )
    public boolean useExtendedCigar = true;

    private SAMFileGATKReadWriter outputWriter;

    @Override
    public boolean requiresReference() { return true; }

    @Override
    public CountingReadFilter makeReadFilter(){
        return new CountingReadFilter(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(OUTPUT, false);
    }

    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext ) {

        final GATKRead sanitizedRead = sanitizeRead(read, referenceContext);

        // Write the read to the output file:
        outputWriter.addRead(read);
    }

    private GATKRead sanitizeRead(final GATKRead read, final ReferenceContext referenceContext) {
        final SimpleInterval readInterval  = new SimpleInterval(read.getContig(), read.getStart(), read.getEnd());
        final byte[]         readBases     = read.getBasesNoCopy();
        final byte[]         readbaseQuals = read.getBaseQualitiesNoCopy();
        final byte[]         refBases      = referenceContext.getBases(readInterval);

        final List<Byte>         newReadBases     = new ArrayList<>();
        final List<Byte>         newBaseQualities = new ArrayList<>();
        final List<CigarElement> newCigar         = new ArrayList<>();

        int readIndex = 0;
        int refIndex = 0;

        // Track our current cigar operator so we can accumulate cigars:
        CigarOperator currentNewCigarOp = null;
        int currentNewCigarOpCount = 0;

        CigarOperator iterCigarOp = null;
        int iterCigarOpCount = 0;

        for ( final CigarElement cigarElement :  read.getCigar().getCigarElements() ) {

            // For these elements we don't have to do anything special:
            if (cigarElement.getOperator().equals(CigarOperator.H) ||
                    cigarElement.getOperator().equals(CigarOperator.N) ||
                    cigarElement.getOperator().equals(CigarOperator.P)) {

                iterCigarOp = cigarElement.getOperator();
                iterCigarOpCount = cigarElement.getLength();
            }
            else if ( cigarElement.getOperator().equals(CigarOperator.S) ) {
                for (int i = 0; i < cigarElement.getLength(); ++i) {
                    newReadBases.add(readBases[readIndex + i]);
                    newBaseQualities.add(readbaseQuals[readIndex + i]);
                }
                iterCigarOp = cigarElement.getOperator();
                iterCigarOpCount = cigarElement.getLength();
            }
            else if ( cigarElement.getOperator().equals(CigarOperator.EQ) ) {
                for (int i = 0; i < cigarElement.getLength(); ++i) {
                    newReadBases.add(readBases[readIndex + i]);
                    newBaseQualities.add(readbaseQuals[readIndex + i]);
                }
                iterCigarOp = cigarElement.getOperator();
                iterCigarOpCount = cigarElement.getLength();
            }
            // For the rest of the elements, we have to do something more:
            else if (cigarElement.getOperator().equals(CigarOperator.M)) {
                for (int i = 0; i < cigarElement.getLength(); ++i) {
                    newReadBases.add(refBases[refIndex + i]);

                    if (readBases[readIndex + i] == refBases[refIndex + i]) {
                        // Since the base is the same as the reference anyway, we keep the qual:
                        newBaseQualities.add(readbaseQuals[readIndex + i]);
                    }
                    else {
                        // Since we replaced the reference base we se teh quality to 60:
                        newBaseQualities.add((byte)refQual);
                    }
                }
                iterCigarOp = useExtendedCigar ? CigarOperator.EQ : CigarOperator.M;
                iterCigarOpCount = cigarElement.getLength();
            }
            else if (cigarElement.getOperator().equals(CigarOperator.X) ||
                    cigarElement.getOperator().equals(CigarOperator.D) ) {

                // For these operators we need to add in the reference bases as matches:
                for (int i = 0; i < cigarElement.getLength(); ++i) {
                    newReadBases.add(refBases[refIndex + i]);
                    // Since we know it's the reference base, we set it to the max quality (60):
                    newBaseQualities.add((byte)refQual);
                }
                iterCigarOp = useExtendedCigar ? CigarOperator.EQ : CigarOperator.M;
                iterCigarOpCount = cigarElement.getLength();
            }
            else if (cigarElement.getOperator().equals(CigarOperator.I)) {
                // Inserted bases are simply removed and ignored.
                iterCigarOp = currentNewCigarOp;
                iterCigarOpCount = 0;
            }
            else {
                throw new UserException.MalformedFile("Unexpected cigar operation: " + cigarElement.toString());
            }

            // Update Cigar:
            if ( iterCigarOp == currentNewCigarOp ) {
                currentNewCigarOpCount += iterCigarOpCount;
            }
            else {
                if ( currentNewCigarOp != null ) {
                    newCigar.add(new CigarElement(currentNewCigarOpCount, currentNewCigarOp));
                }
                currentNewCigarOp = iterCigarOp;
                currentNewCigarOpCount = iterCigarOpCount;
            }

            // Update indices:
            if (cigarElement.getOperator().consumesReferenceBases()) {
                refIndex += cigarElement.getLength();
            }
            if (cigarElement.getOperator().consumesReadBases()) {
                readIndex += cigarElement.getLength();
            }
        }

        // Add in the last cigar element now that we're done iterating:
        if ( iterCigarOp == currentNewCigarOp ) {
            newCigar.add(new CigarElement(currentNewCigarOpCount, currentNewCigarOp));
        }

        // Replace the old read data with the new data:
        read.setCigar(new Cigar(newCigar));

        final Byte[] newBases = newReadBases.toArray(new Byte[newReadBases.size()]);
        read.setBases(ArrayUtils.toPrimitive(newBases));

        final Byte[] newQualitites = newBaseQualities.toArray(new Byte[newBaseQualities.size()]);
        read.setBaseQualities(ArrayUtils.toPrimitive(newQualitites));

        return read;
    }

    @Override
    public void closeTool() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }
}
