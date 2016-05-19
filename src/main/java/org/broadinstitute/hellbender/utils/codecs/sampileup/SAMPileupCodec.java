/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.codecs.sampileup;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.exception.CodecLineParsingException;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.util.ParsingUtils;

import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static org.broadinstitute.gatk.utils.codecs.sampileup.SAMPileupFeature.VariantType;

/**
 * Decoder for SAM pileup data.
 *
 * <p>
 *     From the <a href="http://samtools.sourceforge.net/">SAMTools project documentation</a>:
 * </p>
 * <p>The Pileup format was first used by Tony Cox and Zemin Ning at
 *     the Sanger Institute. It describes the base-pair information at each chromosomal position. This format
 *     facilitates SNP/indel calling and brief alignment viewing by eye. Note that the pileup program has been replaced
 *     in Samtools by mpileup, which produces a slightly different output format by default.
 * </p>

 * <h3>Format</h3>
 * <p>There are two versions of the original pileup format: the current 6-column format produced by Samtools, and the old
 * 10/13-column "consensus" format which could be obtained by using the -c argument, now deprecated. </p>
 * <h4>Simple pileup: 6-column format</h4>
 * <p>
 *     Each line consists of chromosome, 1-based coordinate, reference base, the
 *     number of reads covering the site, read bases and base qualities. At the
 *     read base column, a dot stands for a match to the reference base on the
 *     forward strand, a comma for a match on the reverse strand, `ACGTN' for a mismatch
 *     on the forward strand and `acgtn' for a mismatch on the reverse strand.
 *     A pattern `\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion between
 *     this reference position and the next reference position. The length of the
 *     insertion is given by the integer in the pattern, followed by the inserted sequence.
 * </p>
 * <pre>
 *     seq1 272 T 24  ,.$.....,,.,.,...,,,.,..^+. <<<+;<<<<<<<<<<<=<;<;7<&
 *     seq1 273 T 23  ,.....,,.,.,...,,,.,..A <<<;<<<<<<<<<3<=<<<;<<+
 *     seq1 274 T 23  ,.$....,,.,.,...,,,.,...    7<7;<;<<<<<<<<<=<;<;<<6
 *     seq1 275 A 23  ,$....,,.,.,...,,,.,...^l.  <+;9*<<<<<<<<<=<<:;<<<<
 *     seq1 276 G 22  ...T,,.,.,...,,,.,....  33;+<<7=7<<7<&<<1;<<6<
 *     seq1 277 T 22  ....,,.,.,.C.,,,.,..G.  +7<;<<<<<<<&<=<<:;<<&<
 *     seq1 278 G 23  ....,,.,.,...,,,.,....^k.   %38*<<;<7<<7<=<<<;<<<<<
 *     seq1 279 C 23  A..T,,.,.,...,,,.,..... ;75&<<<<<<<<<=<<<9<<:<<
 * </pre>
 * <p>
 *     See the <a href="http://samtools.sourceforge.net/pileup.shtml">Pileup format documentation</a> for more details.
 * </p>
 *
 * <h4>Consensus pileup: 10/13-column format</h4>
 * <p>The "consensus" or extended pileup consists of the following:
 *  <ul>
 *      <li>original 6 columns as described above</li>
 *      <li>4 extra columns representing consensus values (consensus base, consensus quality, variant quality and maximum mapping quality of the
 * reads covering the sites) for all sites, inserted before the bases and quality strings</li>
 *      <li>3 extra columns indicating counts of reads supporting indels (just for indel sites)</li>
 *  </ul>
 * </p>
 * <h4>Example of consensus pileup for SNP or non-variant sites</h4>
 * <pre>
 *     seq1  60  T  T  66  0  99  13  ...........^~.^~.   9<<55<;<<<<<<
 *     seq1  61  G  G  72  0  99  15  .............^~.^y. (;975&;<<<<<<<<
 *     seq1  62  T  T  72  0  99  15  .$..............    <;;,55;<<<<<<<<
 *     seq1  63  G  G  72  0  99  15  .$.............^~.  4;2;<7:+<<<<<<<
 *     seq1  64  G  G  69  0  99  14  ..............  9+5<;;;<<<<<<<
 *     seq1  65  A  A  69  0  99  14  .$............. <5-2<;;<<<<<<;
 *     seq1  66  C  C  66  0  99  13  .............   &*<;;<<<<<<8<
 *     seq1  67  C  C  69  0  99  14  .............^~.    ,75<.4<<<<<-<<
 *     seq1  68  C  C  69  0  99  14  ..............  576<;7<<<<<8<< *
 * </pre>
 *
 * <h4>Example of consensus pileup for indels</h4>
 * <pre>
 *     Escherichia_coli_K12	3995037	*	*\/*	430	0	37	144	*	+A	143	1	0
 *     Escherichia_coli_K12	3995279	*	*\/*	202	0	36	68	*	+A	67	1	0
 *     Escherichia_coli_K12	3995281	*	*\/*	239	0	36	67	*	-CG	66	1	0
 * </pre>
 * <p>
 *     See <a href="http://samtools.sourceforge.net/cns0.shtml/">Consensus pileup format (deprecated)</a> for more details.
 * </p>
 *
 * <h3>Caveat</h3>
 * <p>Handling of indels is questionable at the moment. Proceed with care.</p>
 *
 *
 * @author Matt Hanna, Geraldine VdAuwera
 * @since 2014
 */
public class SAMPileupCodec extends AsciiFeatureCodec<SAMPileupFeature> {
    // number of tokens expected (6 or 10 are valid, anything else is wrong)
    private static final int basicTokenCount = 6;
    private static final int consensusSNPTokenCount = 10;
    private static final int consensusIndelTokenCount = 13;
    private static final char fldDelim = '\t';
    // allocate once and don't ever bother creating them again:
    private static final String baseA = "A";
    private static final String baseC = "C";
    private static final String baseG = "G";
    private static final String baseT = "T";
    private static final String emptyStr = ""; // we will use this for "reference" allele in insertions

    // codec file extension
    protected static final String FILE_EXT = "samp";

    public SAMPileupCodec() {
        super(SAMPileupFeature.class);
    }

    public SAMPileupFeature decode(String line) {
        //+1 because we want to know if we have more than the max
        String[] tokens = new String[consensusIndelTokenCount+1];

        // split the line
        final int count = ParsingUtils.split(line,tokens,fldDelim);

        SAMPileupFeature feature = new SAMPileupFeature();

        /**
         * Tokens 0, 1, 2 are the same for both formats so they will be interpreted without differentiation.
         * The 10/13-format has 4 tokens inserted after token 2 compared to the 6-format, plus 3 more tokens added at
         * the end for indels. We are currently not making any use of the extra indel tokens.
         *
         * Any token count other than basicTokenCount, consensusSNPTokenCount or consensusIndelTokenCount is wrong.
         */
        final String observedString, bases, quals;

        feature.setChr(tokens[0]);
        feature.setStart(Integer.parseInt(tokens[1]));

        if(tokens[2].length() != 1)  {
            throw new CodecLineParsingException("The SAM pileup line had unexpected base " + tokens[2] + " on line = " + line);
        }
        feature.setRef(tokens[2].charAt(0));

        switch (count) {
            case basicTokenCount:
                bases = tokens[4];
                quals = tokens[5];
                // parsing is pretty straightforward for 6-col format
                if ( feature.getRef() == '*' ) {   // this indicates an indel -- but it shouldn't occur with vanilla 6-col format
                    throw new CodecLineParsingException("Found an indel on line = " + line + " but it shouldn't happen in simple pileup format");
                } else {
                    parseBasesAndQuals(feature, bases, quals);
                    feature.setRefBases(tokens[2].toUpperCase());
                    feature.setEnd(feature.getStart());
                }
                break;
            case consensusSNPTokenCount: // pileup called a SNP or a reference base
                observedString = tokens[3].toUpperCase();
                feature.setFWDAlleles(new ArrayList<String>(2));
                feature.setConsensusConfidence(Double.parseDouble(tokens[4]));
                feature.setVariantConfidence(Double.parseDouble(tokens[5]));
                bases = tokens[8];
                quals = tokens[9];
                // confirm that we have a non-variant, not a mis-parsed indel
                if ( feature.getRef() == '*' ) {
                    throw new CodecLineParsingException("Line parsing of " + line + " says we have a SNP or non-variant but the ref base is '*', which indicates an indel");
                }
                // Parse the SNP or non-variant
                parseBasesAndQuals(feature, bases, quals);
                if ( observedString.length() != 1 ) {
                    throw new CodecLineParsingException( "Line parsing of " + line + " says we have a SNP or non-variant but the genotype token is not a single letter: " + observedString);
                }
                feature.setRefBases(tokens[2].toUpperCase());
                feature.setEnd(feature.getStart());

                char ch = observedString.charAt(0);

                switch ( ch ) {  // record alleles (decompose ambiguous base codes)
                    case 'A': feature.getFWDAlleles().add(baseA); feature.getFWDAlleles().add(baseA); break;
                    case 'C': feature.getFWDAlleles().add(baseC); feature.getFWDAlleles().add(baseC); break;
                    case 'G': feature.getFWDAlleles().add(baseG); feature.getFWDAlleles().add(baseG); break;
                    case 'T': feature.getFWDAlleles().add(baseT); feature.getFWDAlleles().add(baseT); break;
                    case 'M': feature.getFWDAlleles().add(baseA); feature.getFWDAlleles().add(baseC); break;
                    case 'R': feature.getFWDAlleles().add(baseA); feature.getFWDAlleles().add(baseG); break;
                    case 'W': feature.getFWDAlleles().add(baseA); feature.getFWDAlleles().add(baseT); break;
                    case 'S': feature.getFWDAlleles().add(baseC); feature.getFWDAlleles().add(baseG); break;
                    case 'Y': feature.getFWDAlleles().add(baseC); feature.getFWDAlleles().add(baseT); break;
                    case 'K': feature.getFWDAlleles().add(baseG); feature.getFWDAlleles().add(baseT); break;
                }
                if ( feature.getFWDAlleles().get(0).charAt(0) == feature.getRef() && feature.getFWDAlleles().get(1).charAt(0) == feature.getRef() ) feature.setVariantType(VariantType.NONE);
                else {
                    // 	we know that at least one allele is non-ref;
                    // if one is ref and the other is non-ref, or if both are non ref but they are the same (i.e.
                    // homozygous non-ref), we still have 2 allelic variants at the site (e.g. one ref and one nonref)
                    feature.setVariantType(VariantType.SNP);
                    if ( feature.getFWDAlleles().get(0).charAt(0) == feature.getRef() ||
                            feature.getFWDAlleles().get(1).charAt(0) == feature.getRef() ||
                            feature.getFWDAlleles().get(0).equals(feature.getFWDAlleles().get(1))
                            ) feature.setNumNonRef(1);
                    else feature.setNumNonRef(2); // if both observations differ from ref and they are not equal to one another, then we get multiallelic site...
                }
                break;
            case consensusIndelTokenCount:
                observedString = tokens[3].toUpperCase();
                feature.setFWDAlleles(new ArrayList<String>(2));
                feature.setConsensusConfidence(Double.parseDouble(tokens[4]));
                feature.setVariantConfidence(Double.parseDouble(tokens[5]));
                // confirm that we have an indel, not a mis-parsed SNP or non-variant
                if ( feature.getRef() != '*' ) {
                    throw new CodecLineParsingException("Line parsing of " + line + " says we have an indel but the ref base is not '*'");
                }
                // Parse the indel
                parseIndels(observedString,feature) ;
                if ( feature.isDeletion() ) feature.setEnd(feature.getStart()+feature.length()-1);
                else feature.setEnd(feature.getStart()); // if it's not a deletion and we are biallelic, this has got to be an insertion; otherwise the state is inconsistent!!!!
                break;
            default:
                throw new CodecLineParsingException("The SAM pileup line didn't have the expected number of tokens " +
                    "(expected = " + basicTokenCount + " (basic pileup), " + consensusSNPTokenCount +
                    " (consensus pileup for a SNP or non-variant site) or " + consensusIndelTokenCount +
                    " (consensus pileup for an indel); saw = " + count + " on line = " + line + ")");
        }
        return feature;
    }

    /**
     * Can the file be decoded?
     * @param path path the file to test for parsability with this codec
     * @return true if the path has the correct file extension, false otherwise
     */
    @Override
    public boolean canDecode(final String path) { return path.endsWith("." + FILE_EXT); }

    @Override
    public Object readActualHeader(LineIterator lineIterator) {
        // No header for this format
        return null;
    }

    private void parseIndels(String genotype,SAMPileupFeature feature) {
        String [] obs = genotype.split("/"); // get observations, now need to tinker with them a bit

        // if reference allele is among the observed alleles, we will need to take special care of it since we do not have direct access to the reference;
        // if we have an insertion, the "reference" allele is going to be empty; if it it is a deletion, we will deduce the "reference allele" bases
        // from what we have recorded for the deletion allele (e.g. "-CAC")
        boolean hasRefAllele = false;

        for ( int i = 0 ; i < obs.length ; i++ ) {
            if ( obs[i].length() == 1 && obs[i].charAt(0) == '*'  ) {
                hasRefAllele = true;
                feature.getFWDAlleles().add(emptyStr);
                continue;
            }

            String varBases = obs[i].toUpperCase();

            switch ( obs[i].charAt(0) )  {
                case '+':
                    if (!feature.isReference() && !feature.isInsertion()) feature.setVariantType(VariantType.INDEL);
                    else feature.setVariantType(VariantType.INSERTION);
                    feature.setRefBases(emptyStr);
                    break;
                case '-' :
                    if (!feature.isReference() && !feature.isDeletion()) feature.setVariantType(VariantType.INDEL);
                    else feature.setVariantType(VariantType.DELETION);
                    feature.setRefBases(varBases); // remember what was deleted, this will be saved as "reference allele"
                    break;
                default: throw new CodecLineParsingException("Can not interpret observed indel allele record: "+genotype);
            }
            feature.getFWDAlleles().add(varBases);
            feature.setLength(obs[i].length()-1); // inconsistent for non-biallelic indels!!
        }
        if ( hasRefAllele ) {
            // we got at least one ref. allele (out of two recorded)
            if (feature.isReference()) { // both top theories are actually ref allele;
                feature.setNumNonRef(0); // no observations of non-reference allele at all
                feature.setRefBases(emptyStr);
            } else {
                feature.setNumNonRef(1); // hasRefAllele = true, so one allele was definitely ref, hence there is only one left
            }
        } else {
            // we observe two non-ref alleles; they better be the same variant, otherwise the site is not bi-allelic and at the moment we
            // fail to set data in a consistent way.
            if ( feature.getFWDAlleles().get(0).equals(feature.getFWDAlleles().get(1))) feature.setNumNonRef(1);
            else feature.setNumNonRef(2);
        }
        // DONE with indels

    }

    private void parseBasesAndQuals(SAMPileupFeature feature, final String bases, final String quals)
    {
        //System.out.printf("%s%n%s%n", bases, quals);

        // needs to convert the base string with its . and , to the ref base
        StringBuilder baseBuilder = new StringBuilder();
        StringBuilder qualBuilder = new StringBuilder();
        boolean done = false;
        for ( int i = 0, j = 0; i < bases.length() && ! done; i++ ) {
            //System.out.printf("%d %d%n", i, j);
            char c = (char)bases.charAt(i);

            switch ( c ) {
                case '.':   // matches reference
                case ',':   // matches reference
                    baseBuilder.append(feature.getRef());
                    qualBuilder.append(quals.charAt(j++));
                    break;
                case '$':   // end of read
                    break;
                case '*':   // end of indel?
                    j++;
                    break;
                case '^':   // mapping quality
                    i++;
                    break;
                case '+':   // start of indel
                case '-':   // start of indel
                    final Pattern regex = Pattern.compile("([0-9]+).*");             // matches case 1
                    final String rest = bases.substring(i+1);
                    //System.out.printf("sub is %s%n", rest);
                    Matcher match = regex.matcher(rest);
                    if ( ! match.matches() ) {
                        if ( feature.getRef() != '*' )
                            throw new CodecLineParsingException("Bad pileup format: " + bases + " at position " + i);
                        done = true;
                    }
                    else {
                        String g = match.group(1);
                        //System.out.printf("group is %d, match is %s%n", match.groupCount(), g);
                        int l = Integer.parseInt(g);
                        i += l + g.length();    // length of number + that many bases + +/- at the start (included in the next i++)
                        //System.out.printf("remaining is %d => %s%n", l, bases.substring(i+1));
                    }
                    break;
                default:   // non reference base
                    baseBuilder.append(c);
                    qualBuilder.append(quals.charAt(j++));
            }
        }

        feature.setPileupBases(baseBuilder.toString());
        feature.setPileupQuals(qualBuilder.toString());
    }
}
