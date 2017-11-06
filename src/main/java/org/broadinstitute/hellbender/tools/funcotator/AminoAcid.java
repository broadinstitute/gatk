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
    VALINE("Valine","Val","V",new String[]{"GTA","GTC","GTG","GTT"});

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
