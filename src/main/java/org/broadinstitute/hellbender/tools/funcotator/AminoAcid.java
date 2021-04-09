package org.broadinstitute.hellbender.tools.funcotator;

/**
 * Enum to hold the amino acids and their standard codons.
 * Standard codons include the IUPAC base equivalents as found here:
 * https://en.wikipedia.org/wiki/DNA_codon_table
 */
public enum AminoAcid {
    
    ALANINE("Alanine","Ala","A",new String[]{"GCA","GCC","GCG","GCT", "GCN"}),
    ARGANINE("Arganine","Arg","R",new String[]{"AGA","AGG","CGA","CGC","CGG","CGT", "CGN","AGR","CGY","MGR"}),
    ASPARAGINE("Asparagine","Asn","N",new String[]{"AAC","AAT", "AAY"}),
    ASPARTIC_ACID("Aspartic acid","Asp","D",new String[]{"GAT","GAC", "GAY"}),
    CYSTEINE("Cysteine","Cys","C",new String[]{"TGC","TGT", "TGY"}),
    GLUTAMIC_ACID("Glutamic acid","Glu","E",new String[]{"GAA","GAG", "GAR"}),
    GLUTAMINE("Glutamine","Gln","Q",new String[]{"CAA","CAG", "CAR"}),
    GLYCINE("Glycine","Gly","G",new String[]{"GGA","GGC","GGG","GGT", "GGN"}),
    HISTIDINE("Histidine","His","H",new String[]{"CAC","CAT", "CAY"}),
    ISOLEUCINE("Isoleucine","Ile","I",new String[]{"ATA","ATC","ATT", "ATH"}),
    LEUCINE("Leucine","Leu","L",new String[]{"CTA","CTC","CTG","CTT","TTA","TTG", "CTN","CTY","TTR","YTR"}),
    LYSINE("Lysine","Lys","K", new String[]{"AAA","AAG", "AAR"}),
    METHIONINE("Methionine","Met","M",new String[]{"ATG"}),
    PHENYLALANINE("Phenylalanine","Phe","F",new String[]{"TTC","TTT", "TTY"}),
    PROLINE("Proline","Pro","P",new String[]{"CCA","CCC","CCG","CCT", "CCN"}),
    SERINE("Serine","Ser","S",new String[]{"AGC","AGT","TCA","TCC","TCG","TCT", "TCN","AGY"}),
    STOP_CODON("Stop codon","Stop","*",new String[]{"TAA","TAG","TGA", "TRA","TAR"}),
    THREONINE("Threonine","Thr","T",new String[]{"ACA","ACC","ACG","ACT", "ACN"}),
    TRYPTOPHAN("Tryptophan","Trp","W",new String[]{"TGG"}),
    TYROSINE("Tyrosine","Tyr","Y",new String[]{"TAC","TAT", "TAY"}),
    VALINE("Valine","Val","V",new String[]{"GTA","GTC","GTG","GTT", "GTN"}),

    // Need to have an undecodable Amino acid here in case we encounter an IUPAC base that causes a protein
    // sequence to be ambiguous.
    UNDECODABLE("Undecodable Amino Acid", "UNDECODABLE", "?", new String[]{});

    /**
     * The length of a codon in bases.
     */
    public static final int CODON_LENGTH = 3;

    private String[] codons;
    private String fullName;
    private String code;
    private String letter;

    AminoAcid(final String name, final String shortName, final String abbrev, final String[] myCodons) {
        codons = myCodons;
        fullName = name;
        code = shortName;
        letter = abbrev;
    }

    public String getName() {
        return fullName;
    }

    public String getLetter() {
        return letter;
    }

    public String getCode() {
        return code;
    }

    public boolean isStop() {
        return this == STOP_CODON;
    }

    public String toString() {
        return getName();
    }

    public String[] getCodons() {
        return codons;
    }
}
