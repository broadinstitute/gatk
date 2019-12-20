package org.broadinstitute.hellbender.utils.codecs.refseq;

import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.List;

/**
 * the ref seq feature
 */
public class RefSeqFeature implements Transcript, Feature {

    private String transcript_id;
    private int strand;
    private SimpleInterval transcript_interval;
    private SimpleInterval transcript_coding_interval;
    private List<SimpleInterval> exons;
    private String gene_name;
    private List<Integer> exon_frames;
    private String name;

    public RefSeqFeature(SimpleInterval genomeLoc) {
        this.transcript_interval = genomeLoc;
    }

    /** Returns id of the transcript (RefSeq NM_* id) */
    public String getTranscriptId() { return transcript_id; }

    /** Returns coding strand of the transcript, 1 or -1 for positive or negative strand, respectively */
    public int getStrand() { return strand; }

    @Override
    public SimpleInterval getLocation() {
        return transcript_interval;
    }

    /** Returns genomic interval of the coding sequence (does not include UTRs, but still includes introns, since it's a single interval on the DNA) */
    public SimpleInterval getCodingLocation() { return transcript_coding_interval; }

    /** Name of the gene this transcript corresponds to (NOT gene id such as Entrez etc) */
    public String getGeneName() { return gene_name; }

    /** Number of exons in this transcript */
    public int getNumExons() { return exons.size(); }

    /** Genomic location of the n-th exon; throws an exception if n is out of bounds */
    public SimpleInterval getExonLocation(int n) {
        if ( n >= exons.size() || n < 0 ) {
            throw new GATKException("Index out-of-bounds. Transcript has " + exons.size() +" exons; requested: "+n);
        }
        return exons.get(n);
    }

    /** Returns the list of all exons in this transcript, as genomic intervals */
    public List<SimpleInterval> getExons() { return exons; }

//    /** Returns all exons falling ::entirely:: inside an interval **/
//    public List<SimpleInterval> getExonsInInterval( GenomeLoc interval ) {
//        List<SimpleInterval> relevantExons = new ArrayList<>(exons.size());
//        for ( SimpleInterval exon : getExons() ) {
//            if (IntervalUtils.containsP(interval, exon) ) {
//                relevantExons.add(exon);
//            }
//        }
//
//        return relevantExons;
//    }
//
//    /** convenience method; returns the numbers of the exons in the interval **/
//    public List<Integer> getExonNumbersInInterval( Locatable interval ) {
//        List<Integer> numbers = new ArrayList<Integer>();
//        int iNo = 0;
//        for ( Locatable exon : getExons() ) {
//            if ( IntervalUtils.containsP(interval, exon) ) {
//                numbers.add(iNo);
//            }
//            iNo++;
//        }
//
//        return numbers;
//    }

    public String getTranscriptUniqueGeneName() {
        return String.format("%s(%s)",getGeneName(),getTranscriptId());
    }

//    ArrayList<GenomeLoc> exonInRefOrderCache = null;

//    public Integer getSortedOverlapInteger(GenomeLoc position) {
//        int exonNo = -1;
//        ArrayList<GenomeLoc> exonsInReferenceOrder = exonInRefOrderCache != null ? exonInRefOrderCache : new ArrayList<GenomeLoc>(exons);
//        if ( exonInRefOrderCache == null ) {
//            Collections.sort(exonsInReferenceOrder);
//        }
//        exonInRefOrderCache = exonsInReferenceOrder;
//        for ( GenomeLoc exon : exonsInReferenceOrder ) {
//            if ( exon.overlapsP(position) ) {
//                return ++exonNo;
//            }
//            ++exonNo;
//        }
//
//        return -1;
//    }

//    public GenomeLoc getSortedExonLoc(int offset) {
//        ArrayList<GenomeLoc> exonsInReferenceOrder = exonInRefOrderCache != null ? exonInRefOrderCache : new ArrayList<GenomeLoc>(exons);
//        if ( exonInRefOrderCache == null ) {
//            Collections.sort(exonsInReferenceOrder);
//        }
//        exonInRefOrderCache = exonsInReferenceOrder;
//        return exonsInReferenceOrder.get(offset);
//    }

//    /** Returns true if the specified interval 'that' overlaps with the full genomic interval of this transcript */
//    public boolean overlapsP (Locatable that) {
//        return IntervalUtils.overlaps(this, that);
//    }
//
//    /** Returns true if the specified interval 'that' overlaps with the coding genomic interval of this transcript.
//     * NOTE: since "coding interval" is still a single genomic interval, it will not contain UTRs of the outermost exons,
//     * but it will still contain introns and/or exons internal to this genomic locus that are not spliced into this transcript.
//     * @see #overlapsExonP
//     */
//    public boolean overlapsCodingP (Locatable that) {
//        return transcript_coding_interval.overlapsP(that);
//    }
//
    /**
     * Returns a count of the total number of reference bases spanned by gene summary. Will total the length of the exons or
     * if absent, the lengthOnReference for the gene itself.
     *
     * NOTE: This currently makes the assumption that genes do not ever have overlapping exons. I do not know if this is a fair
     *       assumption given extant RefSeqGeneList files.
     */
    public int getTotalExonLength() {
        if (exons.isEmpty()) {
            return getLengthOnReference();
        }
        return exons.stream().mapToInt(Locatable::getLengthOnReference).sum();
    }

    /** Returns true if the specified interval 'that' overlaps with any of the exons actually spliced into this transcript */
    @Override
    public boolean contains(Locatable that) {
        if (exons.isEmpty()) {
            return getLocation().contains(that);
        }
        for ( SimpleInterval exon : exons ) {
            if ( IntervalUtils.overlaps(exon, that) ) {
                return true;
            }
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
            b.append( getLocation().getEnd());  // #5
            b.append('\t');
            b.append( (transcript_coding_interval.getStart() - 1) ); // #6
            b.append('\t');
            b.append( transcript_coding_interval.getEnd());  // #7
            b.append('\t');
            b.append(exons.size()); // #8
            b.append('\t');
            for ( SimpleInterval loc : exons ) { b.append( (loc.getStart()-1) ); b.append(','); } // #9
            b.append('\t');
            for ( SimpleInterval loc : exons ) { b.append( loc.getEnd() ); b.append(','); } // #10
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

    public void setTranscript_interval(SimpleInterval transcript_interval) {
        this.transcript_interval = transcript_interval;
    }

    public void setTranscript_coding_interval(SimpleInterval transcript_coding_interval) {
        this.transcript_coding_interval = transcript_coding_interval;
    }

    public void setExons(List<SimpleInterval> exons) {
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

    @Override
    public String getContig() {
        return transcript_interval.getContig();
    }

    @Override
    public int getStart() {
        return transcript_interval.getStart();
    }

    @Override
    public int getEnd() {
        return transcript_interval.getEnd();
    }
}
