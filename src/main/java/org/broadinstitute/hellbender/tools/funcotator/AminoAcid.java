package org.broadinstitute.hellbender.tools.funcotator;

public enum AminoAcid {
    
    ALANINE("Alanine","Ala","A",new String[]{"GCA","GCC","GCG","GCT"}),
    ARGANINE("Arganine","Arg","R",new String[]{"AGA","AGG","CGA","CGC","CGG","CGT"}),
    ASPARAGINE("Asparagine","Asn","N",new String[]{"AAC","AAT"}),
    ASPARTIC_ACID("Aspartic acid","Asp","D",new String[]{"GAT","GAC"}),
    CYSTEINE("Cysteine","Cys","C",new String[]{"TGC","TGT"}),
    GLUTAMIC_ACID("Glutamic acid","Glu","E",new String[]{"GAA","GAG"}),
    GLUTAMINE("Glutamine","Gln","Q",new String[]{"CAA","CAG"}),
    GLYCINE("Glycine","Gly","G",new String[]{"GGA","GGC","GGG","GGT"}),
    HISTIDINE("Histidine","His","H",new String[]{"CAC","CAT"}),
    ISOLEUCINE("Isoleucine","Ile","I",new String[]{"ATA","ATC","ATT"}),
    LEUCINE("Leucine","Leu","L",new String[]{"CTA","CTC","CTG","CTT","TTA","TTG"}),
    LYSINE("Lysine","Lys","K", new String[]{"AAA","AAG"}),
    METHIONINE("Methionine","Met","M",new String[]{"ATG"}),
    PHENYLALANINE("Phenylalanine","Phe","F",new String[]{"TTC","TTT"}),
    PROLINE("Proline","Pro","P",new String[]{"CCA","CCC","CCG","CCT"}),
    SERINE("Serine","Ser","S",new String[]{"AGC","AGT","TCA","TCC","TCG","TCT"}),
    STOP_CODON("Stop codon","Stop","*",new String[]{"TAA","TAG","TGA"}),
    THREONINE("Threonine","Thr","T",new String[]{"ACA","ACC","ACG","ACT"}),
    TRYPTOPHAN("Tryptophan","Trp","W",new String[]{"TGG"}),
    TYROSINE("Tyrosine","Tyr","Y",new String[]{"TAC","TAT"}),
    VALINE("Valine","Val","V",new String[]{"GTA","GTC","GTG","GTT"}),

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
