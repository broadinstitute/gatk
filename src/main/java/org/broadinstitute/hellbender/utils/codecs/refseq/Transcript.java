/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.hellbender.utils.codecs.refseq;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.HasGenomeLocation;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Sep 22, 2009
 * Time: 5:22:30 PM
 * To change this template use File | Settings | File Templates.
 */
public interface Transcript extends Locatable {

    /** Returns id of the transcript (RefSeq NM_* id) */
    public String getTranscriptId();
    /** Returns coding strand of the transcript, 1 or -1 for positive or negative strand, respectively */
    public int getStrand();
    /** Returns transcript's full genomic interval (includes all exons with UTRs) */
    public SimpleInterval getLocation();
    /** Returns genomic interval of the coding sequence (does not include
     * UTRs, but still includes introns, since it's a single interval on the DNA)
     */
    public SimpleInterval getCodingLocation();
    /** Name of the gene this transcript corresponds to (typically NOT gene id such as Entrez etc,
     * but the implementation can decide otherwise)
     */
    public String getGeneName();
    /** Number of exons in this transcript */
//    public int getNumExons();
//    /** Genomic location of the n-th exon; expected to throw an exception (runtime) if n is out of bounds */
//    public SimpleInterval getExonLocation(int n);

    /** Returns the list of all exons in this transcript, as genomic intervals */
    public List<SimpleInterval> getExons();

//    /** Returns true if the specified interval 'that' overlaps with the full genomic interval of this transcript */
//    public boolean overlapsP(SimpleInterval that);
//
//    /** Returns true if the specified interval 'that' overlaps with the coding genomic interval of this transcript.
//      * NOTE: since "coding interval" is still a single genomic interval, it will not contain UTRs of the outermost exons,
//      * but it will still contain introns and/or exons internal to this genomic locus that are not spliced into this transcript.
//      * @see #overlapsExonP
//      */
//    public boolean overlapsCodingP(SimpleInterval that);
//
//    /** Returns true if the specified interval 'that' overlaps with any of the exons actually spliced into this transcript */
//    public boolean overlapsExonP(SimpleInterval that);


}
