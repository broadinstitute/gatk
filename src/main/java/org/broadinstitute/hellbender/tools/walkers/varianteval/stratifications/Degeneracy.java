package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

/**
 * Experimental stratification by the degeneracy of an amino acid, according to VCF annotation.  Not safe
 */
public class Degeneracy extends VariantStratifier {
    private HashMap<String, HashMap<Integer, String>> degeneracies;

    public Degeneracy(VariantEvalEngine engine) {
        super(engine);

        states.add("1-fold");
        states.add("2-fold");
        states.add("3-fold");
        states.add("4-fold");
        states.add("6-fold");
        states.add("all");

        final HashMap<String, String[]> aminoAcids = new HashMap<>();
        aminoAcids.put("Ile",  new String[]{"ATT", "ATC", "ATA"});
        aminoAcids.put("Leu",  new String[]{"CTT", "CTC", "CTA", "CTG", "TTA", "TTG"});
        aminoAcids.put("Val",  new String[]{"GTT", "GTC", "GTA", "GTG"});
        aminoAcids.put("Phe",  new String[]{"TTT", "TTC"});
        aminoAcids.put("Met",  new String[]{"ATG"});
        aminoAcids.put("Cys",  new String[]{"TGT", "TGC"});
        aminoAcids.put("Ala",  new String[]{"GCT", "GCC", "GCA", "GCG"});
        aminoAcids.put("Gly",  new String[]{"GGT", "GGC", "GGA", "GGG"});
        aminoAcids.put("Pro",  new String[]{"CCT", "CCC", "CCA", "CCG"});
        aminoAcids.put("Thr",  new String[]{"ACT", "ACC", "ACA", "ACG"});
        aminoAcids.put("Ser",  new String[]{"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"});
        aminoAcids.put("Tyr",  new String[]{"TAT", "TAC"});
        aminoAcids.put("Trp",  new String[]{"TGG"});
        aminoAcids.put("Glu",  new String[]{"CAA", "CAG"});
        aminoAcids.put("Asn",  new String[]{"AAT", "AAC"});
        aminoAcids.put("His",  new String[]{"CAT", "CAC"});
        aminoAcids.put("Gln",  new String[]{"GAA", "GAG"});
        aminoAcids.put("Asp",  new String[]{"GAT", "GAC"});
        aminoAcids.put("Lys",  new String[]{"AAA", "AAG"});
        aminoAcids.put("Arg",  new String[]{"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"});
        aminoAcids.put("Stop", new String[]{"TAA", "TAG", "TGA"});

        degeneracies = new HashMap<>();

        for (String aminoAcid : aminoAcids.keySet()) {
            final String[] codons = aminoAcids.get(aminoAcid);

            for (int pos = 0; pos < 3; pos++) {
                HashSet<Character> alleles = new HashSet<Character>();

                for (String codon : codons) {
                    alleles.add(codon.charAt(pos));
                }

                String degeneracy;
                switch (alleles.size()) {
                    case 1:  degeneracy = "1-fold"; break;
                    case 2:  degeneracy = "2-fold"; break;
                    case 3:  degeneracy = "3-fold"; break;
                    case 4:  degeneracy = "4-fold"; break;
                    case 6:  degeneracy = "6-fold"; break;
                    default: degeneracy = "1-fold"; break;
                }

                if (!degeneracies.containsKey(aminoAcid)) {
                    degeneracies.put(aminoAcid, new HashMap<Integer, String>());
                }

                degeneracies.get(aminoAcid).put(pos, degeneracy);
            }
        }
    }

    @Override
    public List<Object> getRelevantStates(final VariantEvalContext context, final VariantContext comp, final String compName, final VariantContext eval, final String evalName, final String sampleName, final String familyName) {
        ArrayList<Object> relevantStates = new ArrayList<>();

        relevantStates.add("all");

        if (eval != null && eval.isVariant()) {
            String type = null;
            String aa = null;
            Integer frame = null;

            if (eval.hasAttribute("refseq.functionalClass")) {
                aa = eval.getAttributeAsString("refseq.variantAA", null);
                frame = eval.getAttributeAsInt("refseq.frame", 0);
            } else if (eval.hasAttribute("refseq.functionalClass_1")) {
                int annotationId = 1;
                String key;

                do {
                    key = String.format("refseq.functionalClass_%d", annotationId);

                    String newtype = eval.getAttributeAsString(key, null);

                    if ( newtype != null &&
                            ( type == null ||
                                    ( type.equals("silent") && !newtype.equals("silent") ) ||
                                    ( type.equals("missense") && newtype.equals("nonsense") ) )
                            ) {
                        type = newtype;

                        String aakey = String.format("refseq.variantAA_%d", annotationId);
                        aa = eval.getAttributeAsString(aakey, null);

                        if (aa != null) {
                            String framekey = String.format("refseq.frame_%d", annotationId);

                            if (eval.hasAttribute(framekey)) {
                                frame = eval.getAttributeAsInt(framekey, 0);
                            }
                        }
                    }

                    annotationId++;
                } while (eval.hasAttribute(key));
            }

            if (aa != null && degeneracies.containsKey(aa) && frame != null) {
                relevantStates.add(degeneracies.get(aa).get(frame));
            }
        }

        return relevantStates;
    }
}
