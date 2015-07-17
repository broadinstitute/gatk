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

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;

import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: 7/4/11
 * Time: 1:09 PM
 *
 * Interface for readers creating diffable objects from a file
 */
public interface DiffableReader {
    @Ensures("result != null")
    /**
     * Return the name of this DiffableReader type.  For example, the VCF reader returns 'VCF' and the
     * bam reader 'BAM'
     */
    public String getName();

    @Ensures("result != null")
    @Requires("file != null")
    /**
     * Read up to maxElementsToRead DiffElements from file, and return them.
     */
    public DiffElement readFromFile(File file, int maxElementsToRead);

    /**
     * Return true if the file can be read into DiffElement objects with this reader. This should
     * be uniquely true/false for all readers, as the system will use the first reader that can read the
     * file.  This routine should never throw an exception.  The VCF reader, for example, looks at the
     * first line of the file for the ##format=VCF4.1 header, and the BAM reader for the BAM_MAGIC value
     * @param file
     * @return
     */
    @Requires("file != null")
    public boolean canRead(File file);
}
