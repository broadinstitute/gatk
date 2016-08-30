package org.broadinstitute.hellbender.utils;


import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;


/**
 * A builder class that expands functionality for SA tags. Each SATagBuilder is associated with a {@link GATKRead} at
 * construction and supports various options of editing existing SA tags as well as adding new SA tags. NOTE: in order
 * for updates to the SA tag for a given read to be added to the read, {@link SATagBuilder#setSATag()} must be called!
 *
 * When an SATagBuilder is constructed, it will automatically parse any existing SA tag strings on the read into
 * processable units. In order to get the SA tag units as reads, simply call {@link SATagBuilder#getArtificialReadsBasedOnSATag(SAMFileHeader)}.
 * Existing SA tag units can be cleared using {@link SATagBuilder#clear}.
 *
 * In order to add a supplementary read to a tag, simply call {@link SATagBuilder#addTag(GATKRead)} or {@link SATagBuilder#addTag(SATagBuilder)}
 * and the handler will add a new section to its SA tag at the begining or the end of the SA string depending if the other
 * read is supplementary or not. In order to edit a read in the Builder, a tool author must call {@link SATagBuilder#removeTag(GATKRead)}
 * or {@link SATagBuilder#removeTag(SATagBuilder)} to remove the old tag and then {@link SATagBuilder#addTag(GATKRead)} to
 * add the edited version of the read back in.
 *
 * Additionally, the static method {@link SATagBuilder#setReadsAsSupplemental(GATKRead, List)} can be used to set a group of reads as being
 * supplementary to each other with correct SA tags. If there are existing SA tags on any of the reads they will be
 * preserved in the operation.
 */
public class SATagBuilder {
    private List<SARead> supplementaryReads;
    final private GATKRead read;
    private SARead thisRead;
    final private static String cigarRe = "\\*|([0-9]+[MIDNSHPX=])+";
    public SATagBuilder(GATKRead read) {
        Utils.nonNull(read);

        this.read = read;
        supplementaryReads = new LinkedList<>();
        parseSATag();
    }

    /**
     * Clears all SA tag units from the SATagBuilder as they exist.
     *
     * @return
     */
    public SATagBuilder clear() {
        supplementaryReads.clear();
        return this;
    }

    /**
     * Sets the stored SA tag information into the SA attribute field of the read encapsulated by the SATagBuilder
     */
    public SATagBuilder setSATag() {
        if (supplementaryReads.isEmpty()) {
            return this;
        }
        read.setAttribute("SA", getTag());
        return this;
    }

    // Class that generates the output tag string for this builder
    private String getTag() {
        StringBuilder stringBuilder = new StringBuilder();
        for (SARead sa : supplementaryReads) {
            stringBuilder.append(sa.toString());
        }
        return stringBuilder.toString();
    }

    // Set a given read as the primary read within its group
    private void setPrimary(boolean p) {
        read.setIsSupplementaryAlignment(!p);
    }

    /**
     * Add a read encapsulated by another SATagBuilder to this tag, if the read is set to be primary, then it will be
     * advanced to the front of the list, otherwise it will be placed at the back of the SA string.
     *
     * @param otherTag a SATagBuilder encapsulating the SARead to add
     * @return
     */
    public SATagBuilder addTag(SATagBuilder otherTag) {
        if (otherTag.isPrimary()) {
            supplementaryReads.add(0, otherTag.getThisRead());
        } else {
            supplementaryReads.add(otherTag.getThisRead());
        }
        return this;
    }

    /**
     * Add a read to the SATag for the encapsulated read. Will add the read either to the front or the back of the list depending
     * on whether or not it is set as supplementary.
     *
     * @param otherRead the read to be added the SA string
     * @return
     */
    public SATagBuilder addTag(GATKRead otherRead) {
        if (!otherRead.isSupplementaryAlignment()) {
            supplementaryReads.add(0, new SARead(otherRead));
        } else {
            supplementaryReads.add(new SARead(otherRead));
        }
        return this;
    }

    /**
     * Method that removes all SATag units from the list corresponding to the provided read based on start location
     *
     * @param read read to be removed from the SATagBuilder
     * @return
     */
    public SATagBuilder removeTag(GATKRead read) {
        return removeTag(read.getContig(), read.getStart());
    }

    /**
     * Method that removes all SATag units from the list corresponding to the provided read based on start location
     *
     * @param otherBuilder SATagBuilder encapsulating a read to be removed from this builder
     * @return
     */
    public SATagBuilder removeTag(SATagBuilder otherBuilder) {
        return removeTag(otherBuilder.read.getContig(), otherBuilder.read.getStart());
    }

    /**
     * Method that searches supplementaryReads for every read matching contig and start and removes them from the list.
     *
     * @param contig contig name of the read to be removed from the array
     * @param start start position of the read to be removed
     * @return
     */
    public SATagBuilder removeTag(String contig, int start) {
        String s = Integer.toString(start);
        for (int i = 0; i < supplementaryReads.size(); i++) {
            SARead saRead = supplementaryReads.get(i);
            if (saRead.contig.equals(contig) && saRead.pos.equals(s)){
                supplementaryReads.remove(i);
                i--;
            }
        }
        return this;
    }

    /**
     * Parses the String SATag for a given read and sets supplementaryReads to
     */
    private void parseSATag() {
        if (!read.hasAttribute("SA")){
            return;
        }
        String[] tags = read.getAttributeAsString("SA").split(";");
        supplementaryReads = Arrays.stream(tags).map(r -> new SARead(r)).collect(Collectors.toList());

    }

    private boolean isPrimary() {
        return !read.isSupplementaryAlignment();
    }

    private SARead getThisRead() {
        if (thisRead==null && read!=null) {
            thisRead = new SARead(read);
        }
        return thisRead;
    }

    /**
     * Returns a list of artificial GATKReads corresponding to the reads described by the SA tag in read.
     * Fills as much information as can be inferred from the original read onto the artificial read copies
     * This function is only used for testing (e.g. testGetArtificialReadsBasedOnSATag, testEmptyNMTagSupport)
     *
     * @param header the header to be used in constructing artificial reads
     */
    public List<GATKRead> getArtificialReadsBasedOnSATag(SAMFileHeader header) {
        List<GATKRead> output = new ArrayList<>(supplementaryReads.size());
        GATKRead readCopy = read.copy();
        readCopy.setAttribute("SA", getTag());
        SAMRecord record = readCopy.convertToSAMRecord(header);

        List<SAMRecord> readRecords = SAMUtils.getOtherCanonicalAlignments(record);
        for (SAMRecord artificialRead : readRecords) {
            output.add(new SAMRecordToGATKReadAdapter(artificialRead));
        }

        return output;
    }

    /**
     * Private subclass corresponding to an SA read "unit" with all of the attributes for a given SA tag section
     */
    private static class SARead {
        // each of the fields that could be in an SA tag
        String contig;
        String pos;
        String strand;
        String cigar;
        String mapQ;
        String NM;

        // builder for SARead that parses the existing SA tag
        private SARead(String SATag) {
            String[] values = SATag.split(",", -1);
            if (values.length != 6) {
                throw new GATKException("Could not parse SATag: "+SATag);
            }
            if (!values[1].equals("*") && Integer.parseInt(values[1]) < 0) {
                throw new GATKException("Could not parse POS in SATag: "+SATag);
            }
            if (!values[3].matches(cigarRe)){
                throw new GATKException("Could not parse cigar in SATag: " + SATag);
            }
            if (!values[4].equals("*") && Integer.parseInt(values[4]) < 0) {
                throw new GATKException("Could not parse MapQ in SATag: " + SATag);
            }
            this.contig = values[0];
            this.pos = values[1];
            this.strand = values[2];
            this.cigar = values[3];
            this.mapQ = values[4];
            this.NM = values[5];
        }

        private SARead(GATKRead read) {
            Utils.nonNull(read);

            this.contig = read.getContig();
            this.pos = read.getStart()+"";
            this.strand = (read.isReverseStrand()?"-":"+");
            this.cigar = read.getCigar().toString();
            this.mapQ = read.getMappingQuality()+"";
            this.NM = (read.hasAttribute("NM")? read.getAttributeAsString("NM"):"*");
        }

        public String toString() {
            return String.format("%s,%s,%s,%s,%s,%s;",
                    ((this.contig!=null)?this.contig:"*"),
                    ((this.pos!=null)?this.pos:"0"),
                    (strand.equals("-")?"-":"+"),
                    ((this.cigar!=null)?this.cigar:"*"),
                    ((this.mapQ!=null)?this.mapQ:"255"),
                    ((this.NM!=null)?this.NM:"*"));
        }
    }


    /**
     * Sets a collection of GATKReads as supplemental reads of the primary read,
     * This tool will set the isSupplemental attribute as well as the 'SA:' tag
     * properly according to the Samspec (example: SA:rname,position,strand,Cigar,mapQ,NM;...)
     * Note: this tool will not unset the primary read as supplemental, futhermore it will simply add to any existing
     * SA tags on the given reads.
     *
     * @param primaryRead the primary read in the set of supplemental reads
     * @param supplementalReads an arbitrarily sorted list of reads to be set as
     *                          supplemental to eachother and the primary read
     *
     * @throws UserException
     */
    public static void setReadsAsSupplemental(GATKRead primaryRead, List<GATKRead> supplementalReads) {
        List<SATagBuilder> supplementalTags = new ArrayList<>();

        supplementalTags.add(new SATagBuilder(primaryRead));

        // create the SA tag string for every read in the sequence
        for (GATKRead read : supplementalReads) {
            read.setIsSupplementaryAlignment(true);
            supplementalTags.add(new SATagBuilder(read));
        }

        for (int i = 0; i < supplementalTags.size(); i++) {
            for (int j = 0; j < supplementalTags.size(); j++) {
                if (i != j) {
                    supplementalTags.get(i).addTag(supplementalTags.get(j));
                }
            }
        }

        for (SATagBuilder read : supplementalTags) {
            read.setSATag();
        }
    }
}

