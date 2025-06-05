package org.broadinstitute.hellbender.tools.gvs.common;

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

    private final int index;
    private final String v37ContigName;
    private static Map<String, ChromosomeEnum> ref37 = new HashMap<>();
    private static final Map<String, ChromosomeEnum> ref38 = new HashMap<>();
    private static final Map<Integer, ChromosomeEnum> decodeValues = new HashMap<>();
    private static final Map<Integer, String> customDecodeValues = new HashMap<>();
    private static Map<String, ChromosomeEnum> currentVersion = null;
    private static Map<String, Integer> customContigMap = null;

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

    public static void setCustomContigMap(Map<String, Integer> customMapping) {
        customContigMap = customMapping;
        for (Map.Entry<String, Integer> customContig : customMapping.entrySet()) {
            customDecodeValues.put(customContig.getValue(), customContig.getKey());
        }
    }

    ChromosomeEnum(int index, String v37identifier) {
        this.index = index;
        this.v37ContigName = v37identifier;
    }

    public static ChromosomeEnum valueOfIndex(int index) {
        return decodeValues.get(index);
    }

    public static String stringValueOfIndex(int index) {
        if (!usingCustomMapping()) {
            return decodeValues.get(index).getContigName();
        } else {
            return customDecodeValues.get(index);
        }
    }

    public static boolean usingCustomMapping() {
        return currentVersion != ref37 && currentVersion != ref38;
    }

    public static ChromosomeEnum valueOfContig(String contig) {
        if (currentVersion == null) {
            throw new RuntimeException("must set reference version");
        } else {
            return currentVersion.get(contig);
        }
    }

    public static Integer integerValueOfContig(String contig) {
        if (!usingCustomMapping()) {
            return valueOfContig(contig).index;
        } else {
            // this is a custom contig mapping that we're dealing with.  Look it up directly
            if (customContigMap == null) {
                throw new RuntimeException("Must supply a custom contig mapping for a non-standard reference");
            }
            if (!customContigMap.containsKey(contig)) {
                System.out.println("Contig " + contig + " not found in custom contig map");
            }
            return customContigMap.get(contig);
        }
    }

    public String getContigName() {
        if (currentVersion==ref37) {
            return v37ContigName;
        } else {
            return name();
        }
    }

    public int getIndex() {
        return index;
    }
}
