package org.broadinstitute.hellbender.tools.walkers.analysis;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.LocusWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Histogram;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Optional;
import java.util.stream.LongStream;

/**
 * Created by tsato on 2/23/18.
 */
@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = CoverageAnalysisProgramGroup.class
)
public class CollectBaseQualityHistogram extends LocusWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "")
    private File outputFile = null;

    private PrintStream outputStream = null;

    final private int MAX_BASE_QUALITY = 40;

    final long[] baseQualityHistogram = new long[MAX_BASE_QUALITY + 1];
    final long[] originalBaseQualityHistogram = new long[MAX_BASE_QUALITY + 1];

    @Override
    public void onTraversalStart(){
        try {
            outputStream = outputFile != null ? new PrintStream(outputFile) : System.out;
        }
        catch ( final FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(outputFile, e);
        }
    }

//    @Override
//    public boolean requiresReference() {
//        return true;
//    }

    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        final ReadPileup readPileup = alignmentContext.getBasePileup();
        Iterator<PileupElement> iterator = readPileup.iterator();

        while (iterator.hasNext()){
            final PileupElement pe = iterator.next();
            final int offset = pe.getOffset();
            final GATKRead read = pe.getRead();
            final byte bq = pe.getQual();
            baseQualityHistogram[bq]++;

            final byte[] oqs = ReadUtils.getOriginalBaseQualities(read);
            if (oqs != null){
                final byte oq = oqs[offset];
                originalBaseQualityHistogram[oq]++;
            }
        }
    }

    @Override
    public Object onTraversalSuccess(){
        outputStream.println(String.format("%s\t%s\t%s", "QUALITY", "COUNT_BEFORE", "COUNT_AFTER"));
        for (int i = 0; i < MAX_BASE_QUALITY + 1; i++){
            outputStream.println(String.format("%d\t%d\t%d", i, baseQualityHistogram[i], originalBaseQualityHistogram[i]));
        }

        Utils.validate(LongStream.of(baseQualityHistogram).sum() == LongStream.of(originalBaseQualityHistogram).sum(),
                "number of BQ's and OQ's in the histogram are different");
        return "SUCCESS";
    }

    @Override
    public void closeTool() {
        if ( outputStream != null ) {
            outputStream.close();
        }
    }
}
