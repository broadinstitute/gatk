package org.broadinstitute.hellbender.tools.walkers.qc;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.PeekableIterator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadQueryNameComparator;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.io.File;
import java.io.FileNotFoundException;

// Should be a GATKTool. ReadWalker filters etc...we don't want any bells and whistles
public class ClipReadsForRSEM extends ReadWalker {

    @Argument(fullName = "output", shortName = "O")
    public File outputTable = new File("duplicate_comaprison_counts.csv");

    @Argument(fullName = "out1")
    public File outSamFile1;

    ReadPair currentReadPair;

    @Override
    public void onTraversalStart(){
        final SAMRecord placeholderRead = new SAMRecord(getHeaderForReads());
        placeholderRead.setReadName("placeholder");
        currentReadPair = new ReadPair(new SAMRecordToGATKReadAdapter(placeholderRead));
    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (read.getName().equals(currentReadPair.getQueryName())){
            currentReadPair.add(read);
            return;
        } else {
            if ()
        }
    }

    public void


    @Override
    public Object onTraversalSuccess(){
       return "SUCCESS";
    }
}
