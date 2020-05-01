package org.broadinstitute.hellbender.tools.variantdb;

import java.util.HashMap;
import java.util.Map;

public enum ChromosomeEnum {
    chr1(1, "1"),
    chr2(2, "2"),
    chr3(3, "3"),
    chr4(4, "4"),
    chr5(5, "5"),
    chr6(6, "6"),
    chr7(7, "7"),
    chr8(8, "8"),
    chr9(9, "9"),
    chr10(10, "10"),
    chr11(11, "11"),
    chr12(12, "12"),
    chr13(13, "13"),
    chr14(14, "14"),
    chr15(15, "15"),
    chr16(16, "16"),
    chr17(17, "17"),
    chr18(18, "18"),
    chr19(19, "19"),
    chr20(20, "20"),
    chr21(21, "21"),
    chr22(22, "22"),
    chrX(23, "X"),
    chrY(24, "Y"),
    chrM(25, "MT");

    int index;
    String v37ContigName;
    private static Map<String, ChromosomeEnum> ref37 = new HashMap<>();
    private static final Map<String, ChromosomeEnum> ref38 = new HashMap<>();
    private static final Map<Integer, ChromosomeEnum> decodeValues = new HashMap<>();
    private static Map<String, ChromosomeEnum> currentVersion = null;

    static {
        for (ChromosomeEnum contig : ChromosomeEnum.values()) {
            ref37.put(contig.v37ContigName, contig);
            ref38.put(contig.name(), contig);
            decodeValues.put(contig.index, contig);
        }
    }

    public static void setRefVersion(String refVersion) {
        if (refVersion.equals("37")) {
            currentVersion = ref37;
        } else {
            currentVersion = ref38;
        }
    }

    ChromosomeEnum(int index, String v37identifier) {
        this.index = index;
        this.v37ContigName = v37identifier;
    }

    public static ChromosomeEnum valueOfIndex(int index) {
        return decodeValues.get(index);
    }

    public static ChromosomeEnum valueOfContig(String contig) {
        if (currentVersion == null) {
            throw new RuntimeException("must set reference version");
        } else {
            return currentVersion.get(contig);
        }
    }

    public String getContigName() {
        if (currentVersion==ref37) {
            return v37ContigName;
        } else {
            return name();
        }
    }

}
