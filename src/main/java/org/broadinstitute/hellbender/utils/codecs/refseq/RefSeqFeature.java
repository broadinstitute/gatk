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

package org.broadinstitute.hellbender.utils.codecs.refseq;

import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.GenomeLoc;

import java.util.*;

/**
 * the ref seq feature
 */
public class RefSeqFeature implements Transcript, Feature {

    private String transcript_id;
    private int strand;
    private GenomeLoc transcript_interval;
    private GenomeLoc transcript_coding_interval;
    private List<GenomeLoc> exons;
    private String gene_name;
    private List<Integer> exon_frames;
    private String name;

    public RefSeqFeature(GenomeLoc genomeLoc) {
        this.transcript_interval = genomeLoc;
    }

    /** Returns id of the transcript (RefSeq NM_* id) */
    public String getTranscriptId() { return transcript_id; }

    /** Returns coding strand of the transcript, 1 or -1 for positive or negative strand, respectively */
    public int getStrand() { return strand; }

    /** Returns transcript's full genomic interval (includes all exons with UTRs) */
    public GenomeLoc getLocation() {
        return transcript_interval;
    }

    /** Returns genomic interval of the coding sequence (does not include UTRs, but still includes introns, since it's a single interval on the DNA) */
    public GenomeLoc getCodingLocation() { return transcript_coding_interval; }

    /** Name of the gene this transcript corresponds to (NOT gene id such as Entrez etc) */
    public String getGeneName() { return gene_name; }

    /** Number of exons in this transcript */
    public int getNumExons() { return exons.size(); }

    /** Genomic location of the n-th exon (0-based); throws an exception if n is out of bounds */
    public GenomeLoc getExonLocation(int n) {
        if ( n >= exons.size() || n < 0 ) throw new GATKException("Index out-of-bounds. Transcript has " + exons.size() +" exons; requested: "+n);
        return exons.get(n);
    }

    /** Returns the list of all exons in this transcript, as genomic intervals */
    public List<GenomeLoc> getExons() { return exons; }

    /** Returns all exons falling ::entirely:: inside an interval **/
    public List<GenomeLoc> getExonsInInterval( GenomeLoc interval ) {
        List<GenomeLoc> relevantExons = new ArrayList<GenomeLoc>(exons.size());
        for ( GenomeLoc exon : getExons() ) {
            if ( interval.containsP(exon) ) {
                relevantExons.add(exon);
            }
        }

        return relevantExons;
    }

    /** convenience method; returns the numbers of the exons in the interval (0-based) **/
    public List<Integer> getExonNumbersInInterval( GenomeLoc interval ) {
        List<Integer> numbers = new ArrayList<Integer>();
        int iNo = 0;
        for ( GenomeLoc exon : getExons() ) {
            if ( interval.containsP(exon) ) {
                numbers.add(iNo);
            }
            iNo++;
        }

        return numbers;
    }

    public String getTranscriptUniqueGeneName() {
        return String.format("%s(%s)",getGeneName(),getTranscriptId());
    }

    public String getOverlapString(GenomeLoc position) {
        boolean is_exon = false;
        StringBuilder overlapString = new StringBuilder();
        int exonNo = 1;

        for ( GenomeLoc exon : exons ) {
            if ( exon.containsP(position) ) {
                overlapString.append(String.format("exon_%d",exonNo));
                is_exon = true;
                break;
            }
            exonNo ++;
        }

        if ( ! is_exon ) {
            if ( overlapsCodingP(position) ) {
                overlapString.append("Intron");
            } else {
                overlapString.append("UTR");
            }
        }

        return overlapString.toString();
    }

    ArrayList<GenomeLoc> exonInRefOrderCache = null;

    public Integer getSortedOverlapInteger(GenomeLoc position) {
        int exonNo = -1;
        ArrayList<GenomeLoc> exonsInReferenceOrder = exonInRefOrderCache != null ? exonInRefOrderCache : new ArrayList<GenomeLoc>(exons);
        if ( exonInRefOrderCache == null ) {
            Collections.sort(exonsInReferenceOrder);
        }
        exonInRefOrderCache = exonsInReferenceOrder;
        for ( GenomeLoc exon : exonsInReferenceOrder ) {
            if ( exon.overlapsP(position) ) {
                return ++exonNo;
            }
            ++exonNo;
        }

        return -1;
    }

    public GenomeLoc getSortedExonLoc(int offset) {
        ArrayList<GenomeLoc> exonsInReferenceOrder = exonInRefOrderCache != null ? exonInRefOrderCache : new ArrayList<GenomeLoc>(exons);
        if ( exonInRefOrderCache == null ) {
            Collections.sort(exonsInReferenceOrder);
        }
        exonInRefOrderCache = exonsInReferenceOrder;
        return exonsInReferenceOrder.get(offset);
    }

    /** Returns true if the specified interval 'that' overlaps with the full genomic interval of this transcript */
    public boolean overlapsP (GenomeLoc that) {
        return getLocation().overlapsP(that);
    }

    /** Returns true if the specified interval 'that' overlaps with the coding genomic interval of this transcript.
     * NOTE: since "coding interval" is still a single genomic interval, it will not contain UTRs of the outermost exons,
     * but it will still contain introns and/or exons internal to this genomic locus that are not spliced into this transcript.
     * @see #overlapsExonP
     */
    public boolean overlapsCodingP (GenomeLoc that) {
        return transcript_coding_interval.overlapsP(that);
    }

    /** Returns true if the specified interval 'that' overlaps with any of the exons actually spliced into this transcript */
    public boolean overlapsExonP (GenomeLoc that) {
        for ( GenomeLoc e : exons ) {
            if ( e.overlapsP(that) ) return true;
        }
        return false;
    }

    public String toString() {
        StringBuilder b = new StringBuilder("000\t"); // first field is unused but required in th ecurrent format; just set to something
        b.append(transcript_id);   // #1
        b.append('\t');
        b.append(getLocation().getContig()); // #2
        b.append('\t');
        b.append( (strand==1?'+':'-') ); // #3
        b.append('\t');
        b.append( (getLocation().getStart() - 1) ); // #4
        b.append('\t');
        b.append( getLocation().getStop());  // #5
        b.append('\t');
        b.append( (transcript_coding_interval.getStart() - 1) ); // #6
        b.append('\t');
        b.append( transcript_coding_interval.getStop());  // #7
        b.append('\t');
        b.append(exons.size()); // #8
        b.append('\t');
        for ( GenomeLoc loc : exons ) { b.append( (loc.getStart()-1) ); b.append(','); } // #9
        b.append('\t');
        for ( GenomeLoc loc : exons ) { b.append( loc.getStop() ); b.append(','); } // #10
        b.append("\t0\t"); // # 11 - unused?
        b.append(gene_name); // # 12
        b.append("\tcmpl\tcmpl\t"); // #13, #14 - unused?
        for ( Integer f : exon_frames ) { b.append( f ); b.append(','); } // #15

        return b.toString();
    }

    public void setTranscript_id(String transcript_id) {
        this.transcript_id = transcript_id;
    }

    public void setStrand(int strand) {
        this.strand = strand;
    }

    public void setTranscript_interval(GenomeLoc transcript_interval) {
        this.transcript_interval = transcript_interval;
    }

    public void setTranscript_coding_interval(GenomeLoc transcript_coding_interval) {
        this.transcript_coding_interval = transcript_coding_interval;
    }

    public void setExons(List<GenomeLoc> exons) {
        this.exons = exons;
    }

    public void setGene_name(String gene_name) {
        this.gene_name = gene_name;
    }

    public void setExon_frames(List<Integer> exon_frames) {
        this.exon_frames = exon_frames;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getChr() {
        return transcript_interval.getContig();
    }

    public int getStart() {
        return transcript_interval.getStart();
    }

    public int getEnd() {
        return transcript_interval.getStop();
    }
}
