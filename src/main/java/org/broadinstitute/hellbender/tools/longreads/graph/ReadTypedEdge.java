package org.broadinstitute.hellbender.tools.longreads.graph;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseEdge;

import java.util.regex.Pattern;

public class ReadTypedEdge extends BaseEdge {

    private static final long serialVersionUID = 0x1337L;

    private static final Pattern CCS_READ_PATTERN       = Pattern.compile("/ccs[ \t]*$");
    private static final Pattern RECLAIMED_READ_PATTERN = Pattern.compile("^m[0-9]*.*/[0-9]*/[0-9]*_[0-9]*_[0-9]*_[0-9]*_[0-9]*[ \t]*$");

    private static final String CCS_LABEL = "CCS";
    private static final String RECLAIMED_LABEL = "RECLAIMED";
    private static final String RAW_READ_LABEL = "RAW";

    private String readType;

    public ReadTypedEdge(final String readType) {
        super(false, 1);
        this.readType = readType;
    }

    public String getReadType() {
        return readType;
    }

    public void setReadType(final String readType) {
        this.readType = readType;
    }

    /**
     * Creates a read label string based on the given read name.
     * @param vertex {@link AlignedBaseVertex} containing the name of the read to label.
     * @return A {@link String} containing a read label.
     */
    public static String getReadLabel(final AlignedBaseVertex vertex) {
        return getReadLabel(vertex.getReadName());
    }

    /**
     * Creates a read label string based on the given read name.
     * @param readName {@link String} containing the name of the read to label.
     * @return A {@link String} containing a read label.
     */
    public static String getReadLabel(final String readName) {
//            // TODO: DEBUGGING TO BE REMOVED AT CODE REVIEW TIME
//            System.out.println("ReadName = " + readName);
//            System.out.println("CCS_PATTERN_FIND = " + CCS_READ_PATTERN.matcher(readName).find());
//            System.out.println("RECLAIMED_READ_PATTERN_FIND = " + RECLAIMED_READ_PATTERN.matcher(readName).find());

        if ( CCS_READ_PATTERN.matcher(readName).find() ) {
            // This is a CCS read!
            return CCS_LABEL;
        }
        else if ( RECLAIMED_READ_PATTERN.matcher(readName).find() ) {
            // This is a RECLAIMED read!
            return RECLAIMED_LABEL;
        }
        else {
            // This is a NORMAL read!
            return RAW_READ_LABEL;
        }
    }

    @Override
    public String getDotExtraInfo() {
        return "readType=\"" + readType + "\"";
    }

    @Override
    public BaseEdge copy() {
        return new ReadTypedEdge( readType );
    }

    @Override
    public String getGexfAttributesString() {
        return "<attvalue for=\"0\" value=\"" + readType + "\" />";
    }
}
