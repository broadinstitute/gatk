/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.hellbender.utils.diffengine;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.BlockCompressedInputStream;

import java.io.*;
import java.util.Arrays;


/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: 7/4/11
 * Time: 1:09 PM
 *
 * Class implementing diffnode reader for VCF
 */
public class BAMDiffableReader implements DiffableReader {
    @Override
    public String getName() { return "BAM"; }

    @Override
    public DiffElement readFromFile(File file, int maxElementsToRead) {
        final SAMFileReader reader = new SAMFileReader(file, null); // null because we don't want it to look for the index
        reader.setValidationStringency(ValidationStringency.SILENT);

        DiffNode root = DiffNode.rooted(file.getName());
        SAMRecordIterator iterator = reader.iterator();

        int count = 0;
        while ( iterator.hasNext() ) {
            final SAMRecord record = iterator.next();

            // name is the read name + first of pair
            String name = record.getReadName().replace('.', '_');
            if ( record.getReadPairedFlag() ) {
                name += record.getFirstOfPairFlag() ? "_1" : "_2";
            }

            DiffNode readRoot = DiffNode.empty(name, root);

            // add fields
            readRoot.add("NAME", record.getReadName());
            readRoot.add("FLAGS", record.getFlags());
            readRoot.add("RNAME", record.getReferenceName());
            readRoot.add("POS", record.getAlignmentStart());
            readRoot.add("MAPQ", record.getMappingQuality());
            readRoot.add("CIGAR", record.getCigarString());
            readRoot.add("RNEXT", record.getMateReferenceName());
            readRoot.add("PNEXT", record.getMateAlignmentStart());
            readRoot.add("TLEN", record.getInferredInsertSize());
            readRoot.add("SEQ", record.getReadString());
            readRoot.add("QUAL", record.getBaseQualityString());

            for ( SAMRecord.SAMTagAndValue xt : record.getAttributes() ) {
                readRoot.add(xt.tag, xt.value);
            }

            // add record to root
            if ( ! root.hasElement(name) )
                // protect ourselves from malformed files
                root.add(readRoot);
            count += readRoot.size();
            if ( count > maxElementsToRead && maxElementsToRead != -1)
                break;
        }

        reader.close();

        return root.getBinding();
    }

    @Override
    public boolean canRead(File file) {
        final byte[] BAM_MAGIC = "BAM\1".getBytes();
        final byte[] buffer = new byte[BAM_MAGIC.length];
        try {
            InputStream fstream = new BufferedInputStream(new FileInputStream(file));
            if ( !BlockCompressedInputStream.isValidFile(fstream) )
                return false;
            final BlockCompressedInputStream BCIS = new BlockCompressedInputStream(fstream);
            BCIS.read(buffer, 0, BAM_MAGIC.length);
            BCIS.close();
            return Arrays.equals(buffer, BAM_MAGIC);
        } catch ( IOException e ) {
            return false;
        } catch ( htsjdk.samtools.FileTruncatedException e ) {
            return false;
        }
    }
}
