/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.sam.SamComparison;

import java.io.File;
import java.util.*;

/**
 * Wrapper CLP for SamComparison.
 */
@CommandLineProgramProperties(
        usage = "USAGE: CompareSAMs <SAMFile1> <SAMFile2>\n" +
                "Compares the headers of the two input SAM or BAM files, and, if possible, the SAMRecords. " +
                "For SAMRecords, compares only the readUnmapped flag, reference name, start position and strand. " +
                "Reports the number of SAMRecords that match, differ in alignment, are mapped in only one input, " +
                "or are missing in one of the files",
        usageShort = "Compares two input SAM or BAM files",
        programGroup = ReadProgramGroup.class
)
public class CompareSAMs extends PicardCommandLineProgram {

    @PositionalArguments(minElements = 2, maxElements = 2)
    public List<File> samFiles;

    @Override
    protected Object doWork() {
        SamReaderFactory factory = SamReaderFactory.makeDefault();
        SamReader sam1 = factory.referenceSequence(REFERENCE_SEQUENCE).open(samFiles.get(0));
        SamReader sam2 = factory.referenceSequence(REFERENCE_SEQUENCE).open(samFiles.get(1));
        SamComparison comparison = new SamComparison(sam1, sam2);
        comparison.printReport();
        if (comparison.areEqual()) {
            System.out.println("SAM files match.");
        } else {
            System.out.println("SAM files differ.");
        }
        CloserUtil.close(sam1);
        CloserUtil.close(sam2);
        return null;
    }
}
