package org.broadinstitute.hellbender.utils.variant;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.TribbleException;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFSimpleHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.walkers.genotyper.*;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.io.Serializable;
import java.util.*;
import java.util.function.BiFunction;
import java.util.stream.Collectors;

public final class GATKVariantContextUtils {

    private static final Logger logger = LogManager.getLogger(GATKVariantContextUtils.class);

    public static final String MERGE_FILTER_PREFIX = "filterIn";
    public static final String MERGE_REF_IN_ALL = "ReferenceInAlln";
    public static final String MERGE_FILTER_IN_ALL = "FilteredInAll";
    public static final String MERGE_INTERSECTION = "Intersection";

    public static final int DEFAULT_PLOIDY = HomoSapiensConstants.DEFAULT_PLOIDY;

    private static final GenotypeLikelihoodCalculators GL_CALCS = new GenotypeLikelihoodCalculators();

    public static final double SUM_GL_THRESH_NOCALL = -0.1; // if sum(gl) is bigger than this threshold, we treat GL's as non-informative and will force a no-call.

    public static boolean isInformative(final double[] gls) {
        return MathUtils.sum(gls) < GATKVariantContextUtils.SUM_GL_THRESH_NOCALL;
    }

    /**
     * A method that constructs a mutable set of VCF header lines containing the tool name, version, date and command line,
     * and return to the requesting tool.
     * @param toolkitShortName  short name for the tool kit, e.g. "gatk"
     * @param toolName          name of the tool
     * @param versionString     version of the tool
     * @param dateTime          date and time at invocation of the tool
     * @param cmdLine           command line (e.g. options, flags) when the tool is invoked.
     * @return A mutable set of VCF header lines containing the tool name, version, date and command line.
     */
    public static Set<VCFHeaderLine> getDefaultVCFHeaderLines(final String toolkitShortName, final String toolName,
                                                              final String versionString, final String dateTime,
                                                              final String cmdLine) {
        final Set<VCFHeaderLine> defaultVCFHeaderLines = new HashSet<>();
        final Map<String, String> simpleHeaderLineMap = new HashMap<>(4);
        simpleHeaderLineMap.put("ID", toolName);
        simpleHeaderLineMap.put("Version", versionString);
        simpleHeaderLineMap.put("Date", dateTime);
        simpleHeaderLineMap.put("CommandLine", cmdLine);
        defaultVCFHeaderLines.add(new VCFHeaderLine("source", toolName));
        defaultVCFHeaderLines.add(new VCFSimpleHeaderLine(String.format("%sCommandLine", toolkitShortName), simpleHeaderLineMap));
        return defaultVCFHeaderLines;
    }

    /**
     * Creates a VariantContextWriter whose outputFile type is based on the extension of the output file name.
     * The default options set by VariantContextWriter are cleared before applying ALLOW_MISSING_FIELDS_IN_HEADER (if
     * <code>lenientProcessing</code> is set), followed by the set of options specified by any <code>options</code> args.
     *
     * @param outFile output File for this writer. May not be null.
     * @param referenceDictionary required if on the fly indexing is set, otherwise can be null
     * @param createMD5 true if an md5 file should be created
     * @param options variable length list of additional Options to be set for this writer
     * @returns VariantContextWriter must be closed by the caller
     */
    public static VariantContextWriter createVCFWriter(
            final File outFile,
            final SAMSequenceDictionary referenceDictionary,
            final boolean createMD5,
            final Options... options)
    {
        Utils.nonNull(outFile);

        VariantContextWriterBuilder vcWriterBuilder =
                new VariantContextWriterBuilder().clearOptions().setOutputFile(outFile);

        if (VariantContextWriterBuilder.OutputType.UNSPECIFIED == getVariantFileTypeFromExtension(outFile)) {
            // the only way the user has to specify an output type is by file extension, and htsjdk
            // throws if it can't map the file extension to a known vcf type, so fallback to a default
            // of VCF
            logger.warn(String.format(
                    "Can't determine output variant file format from output file extension \"%s\". Defaulting to VCF.",
                    FilenameUtils.getExtension(outFile.getPath())));
            vcWriterBuilder = vcWriterBuilder.setOutputFileType(VariantContextWriterBuilder.OutputType.VCF);
        }

        if (createMD5) {
            vcWriterBuilder.setCreateMD5();
        }

        if (null != referenceDictionary) {
            vcWriterBuilder = vcWriterBuilder.setReferenceDictionary(referenceDictionary);
        }

        for (Options opt : options) {
            vcWriterBuilder = vcWriterBuilder.setOption(opt);
        }

        return vcWriterBuilder.build();
    }

    // Determine the variant file type from the file extension. Htsjdk has similar code, when
    // https://github.com/broadinstitute/gatk/issues/2128 is fixed we should eliminate this code
    // and use the htsjdk method.
    private static VariantContextWriterBuilder.OutputType getVariantFileTypeFromExtension(final File outputFile) {
        final String extension = FilenameUtils.getExtension(outputFile.getPath()).toLowerCase();

        if (extension.equals(VcfUtils.VCF_FILE_EXTENSION)) {
            return VariantContextWriterBuilder.OutputType.VCF;
        } else if (extension.equals(VcfUtils.BCF_FILE_EXTENSION)) {
            return VariantContextWriterBuilder.OutputType.BCF;
        } else if (IOUtil.hasBlockCompressedExtension(outputFile.getPath())) {
            return VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF;
        }
        return VariantContextWriterBuilder.OutputType.UNSPECIFIED;
    }

    /**
     * Diploid NO_CALL allele list...
     *
     * @deprecated you should use {@link #noCallAlleles(int)} instead. It indicates the presence of a hardcoded diploid assumption which is bad.
     */
    @Deprecated
    public final static List<Allele> NO_CALL_ALLELES = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);

    private static boolean hasPLIncompatibleAlleles(final Collection<Allele> alleleSet1, final Collection<Allele> alleleSet2) {
        final Iterator<Allele> it1 = alleleSet1.iterator();
        final Iterator<Allele> it2 = alleleSet2.iterator();

        while ( it1.hasNext() && it2.hasNext() ) {
            final Allele a1 = it1.next();
            final Allele a2 = it2.next();
            if ( ! a1.equals(a2) )
                return true;
        }

        // by this point, at least one of the iterators is empty.  All of the elements
        // we've compared are equal up until this point.  But it's possible that the
        // sets aren't the same size, which is indicated by the test below.  If they
        // are of the same size, though, the sets are compatible
        return it1.hasNext() || it2.hasNext();
    }

    /**
     * Determines the common reference allele
     *
     * @param VCs    the list of VariantContexts
     * @param loc    if not null, ignore records that do not begin at this start location
     * @return possibly null Allele
     */
    public static Allele determineReferenceAllele(final List<VariantContext> VCs, final Locatable loc) {
        Allele ref = null;

        for ( final VariantContext vc : VCs ) {
            if ( contextMatchesLoc(vc, loc) ) {
                final Allele myRef = vc.getReference();
                if ( ref == null || ref.length() < myRef.length() )
                    ref = myRef;
                else if ( ref.length() == myRef.length() && ! ref.equals(myRef) )
                    throw new TribbleException(String.format("The provided variant file(s) have inconsistent references for the same position(s) at %s:%d, %s vs. %s", vc.getContig(), vc.getStart(), ref, myRef));
            }
        }

        return ref;
    }

    /**
     * Calculates the total ploidy of a variant context as the sum of all plodies across genotypes.
     * @param vc the target variant context.
     * @param defaultPloidy the default ploidy to be assume when there is no ploidy information for a genotype.
     * @return never {@code null}.
     */
    public static int totalPloidy(final VariantContext vc, final int defaultPloidy) {
        Utils.nonNull(vc, "the vc provided cannot be null");
        Utils.validateArg(defaultPloidy >= 0, "the default ploidy must 0 or greater");
        return vc.getGenotypes().stream().mapToInt(Genotype::getPloidy)
                .map(p -> p <= 0 ? defaultPloidy : p).sum();
    }

    /**
     * Returns the type of the substitution (transition vs transversion) represented by this variant. Only applicable to bi allelic SNPs.
     */
    public static BaseUtils.BaseSubstitutionType getSNPSubstitutionType(final VariantContext context) {
        Utils.nonNull(context);
        if (!context.isSNP() || !context.isBiallelic()) {
            throw new IllegalArgumentException("Requested SNP substitution type for bialleic non-SNP " + context);
        }
        return BaseUtils.SNPSubstitutionType(context.getReference().getBases()[0], context.getAlternateAllele(0).getBases()[0]);
    }

    /**
     * If this is a BiAllelic SNP, is it a transition?
     */
    public static boolean isTransition(final VariantContext context) {
        Utils.nonNull(context);
        return getSNPSubstitutionType(context) == BaseUtils.BaseSubstitutionType.TRANSITION;
    }


    /**
     * Returns a homozygous call allele list given the only allele and the ploidy.
     *
     * @param allele the only allele in the allele list.
     * @param ploidy the ploidy of the resulting allele list.
     *
     * @throws IllegalArgumentException if {@code allele} is {@code null} or ploidy is negative.
     *
     * @return never {@code null}.
     */
    public static List<Allele> homozygousAlleleList(final Allele allele, final int ploidy) {
        Utils.nonNull(allele);
        Utils.validateArg(ploidy >= 0, "ploidy cannot be negative");

        // Use a tailored inner class to implement the list:
        return Collections.nCopies(ploidy,allele);
    }

    /**
     * Add the genotype call (GT) field to GenotypeBuilder using the requested {@link GenotypeAssignmentMethod}
     *
     * @param gb the builder where we should put our newly called alleles, cannot be null
     * @param assignmentMethod the method to use to do the assignment, cannot be null
     * @param genotypeLikelihoods a vector of likelihoods to use if the method requires PLs, should be log10 likelihoods, cannot be null
     * @param allelesToUse the alleles with respect to which the likelihoods are defined
     */
    public static void makeGenotypeCall(final int ploidy,
                                        final GenotypeBuilder gb,
                                        final GenotypeAssignmentMethod assignmentMethod,
                                        final double[] genotypeLikelihoods,
                                        final List<Allele> allelesToUse,
                                        final List<Allele> originalGT) {
        if(originalGT == null && assignmentMethod == GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL) {
            throw new IllegalArgumentException("origianlGT cannot be null if assignmentMethod is BEST_MATCH_TO_ORIGINAL");
        }
        if (assignmentMethod == GenotypeAssignmentMethod.SET_TO_NO_CALL) {
            gb.alleles(noCallAlleles(ploidy)).noGQ();
        } else if (assignmentMethod == GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN) {
            if ( genotypeLikelihoods == null || !isInformative(genotypeLikelihoods) ) {
                gb.alleles(noCallAlleles(ploidy)).noGQ();
            } else {
                final int maxLikelihoodIndex = MathUtils.maxElementIndex(genotypeLikelihoods);
                final GenotypeLikelihoodCalculator glCalc = GL_CALCS.getInstance(ploidy, allelesToUse.size());
                final GenotypeAlleleCounts alleleCounts = glCalc.genotypeAlleleCountsAt(maxLikelihoodIndex);

                gb.alleles(alleleCounts.asAlleleList(allelesToUse));
                final int numAltAlleles = allelesToUse.size() - 1;
                if ( numAltAlleles > 0 ) {
                    gb.log10PError(GenotypeLikelihoods.getGQLog10FromLikelihoods(maxLikelihoodIndex, genotypeLikelihoods));
                }
            }
        } else if (assignmentMethod == GenotypeAssignmentMethod.SET_TO_NO_CALL_NO_ANNOTATIONS) {
            gb.alleles(noCallAlleles(ploidy)).noGQ().noAD().noPL().noAttributes();
        } else if (assignmentMethod == GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL) {
            final List<Allele> best = new LinkedList<>();
            final Allele ref = allelesToUse.get(0);
            for (final Allele originalAllele : originalGT) {
                best.add((allelesToUse.contains(originalAllele) || originalAllele.isNoCall()) ? originalAllele : ref);
            }
            gb.alleles(best);
        }
    }

    public static void makeGenotypeCall(final int ploidy,
                                        final GenotypeBuilder gb,
                                        final GenotypeAssignmentMethod assignmentMethod,
                                        final double[] genotypeLikelihoods,
                                        final List<Allele> allelesToUse){
        makeGenotypeCall(ploidy,gb,assignmentMethod,genotypeLikelihoods,allelesToUse,null);
    }

    /**
     * Return the rightmost variant context in maybeOverlapping that overlaps curPos
     *
     * @param curPos non-null genome loc
     * @param maybeOverlapping a collection of variant contexts that might overlap curPos
     * @return the rightmost VariantContext, or null if none overlaps
     */
    public static VariantContext getOverlappingVariantContext(final Locatable curPos, final Collection<VariantContext> maybeOverlapping) {
        VariantContext overlaps = null;
        for ( final VariantContext vc : maybeOverlapping ) {
            if ( curPos.overlaps(vc) ) {
                if ( overlaps == null || vc.getStart() > overlaps.getStart() ) {
                    overlaps = vc;
                }
            }
        }
        return overlaps;
    }

    public enum GenotypeMergeType {
        /**
         * Make all sample genotypes unique by file. Each sample shared across RODs gets named sample.ROD.
         */
        UNIQUIFY,
        /**
         * Take genotypes in priority order (see the priority argument).
         */
        PRIORITIZE,
        /**
         * Take the genotypes in any order.
         */
        UNSORTED,
        /**
         * Require that all samples/genotypes be unique between all inputs.
         */
        REQUIRE_UNIQUE
    }

    public enum FilteredRecordMergeType {
        /**
         * Union - leaves the record if any record is unfiltered.
         */
        KEEP_IF_ANY_UNFILTERED,
        /**
         * Requires all records present at site to be unfiltered. VCF files that don't contain the record don't influence this.
         */
        KEEP_IF_ALL_UNFILTERED,
        /**
         * If any record is present at this site (regardless of possibly being filtered), then all such records are kept and the filters are reset.
         */
        KEEP_UNCONDITIONAL
    }

    /**
     * Returns true iff VC is an non-complex indel where every allele represents an expansion or
     * contraction of a series of identical bases in the reference.
     *
     * For example, suppose the ref bases are CTCTCTGA, which includes a 3x repeat of CTCTCT
     *
     * If VC = -/CT, then this function returns true because the CT insertion matches exactly the
     * upcoming reference.
     * If VC = -/CTA then this function returns false because the CTA isn't a perfect match
     *
     * Now consider deletions:
     *
     * If VC = CT/- then again the same logic applies and this returns true
     * The case of CTA/- makes no sense because it doesn't actually match the reference bases.
     *
     * The logic of this function is pretty simple.  Take all of the non-null alleles in VC.  For
     * each insertion allele of n bases, check if that allele matches the next n reference bases.
     * For each deletion allele of n bases, check if this matches the reference bases at n - 2 n,
     * as it must necessarily match the first n bases.  If this test returns true for all
     * alleles you are a tandem repeat, otherwise you are not.
     *
     * @param vc
     * @param refBasesStartingAtVCWithPad not this is assumed to include the PADDED reference
     * @return
     */
    public static boolean isTandemRepeat(final VariantContext vc, final byte[] refBasesStartingAtVCWithPad) {
        final String refBasesStartingAtVCWithoutPad = new String(refBasesStartingAtVCWithPad).substring(1);
        if ( ! vc.isIndel() ) // only indels are tandem repeats
            return false;

        final Allele ref = vc.getReference();

        for ( final Allele allele : vc.getAlternateAlleles() ) {
            if ( ! isRepeatAllele(ref, allele, refBasesStartingAtVCWithoutPad) )
                return false;
        }

        // we've passed all of the tests, so we are a repeat
        return true;
    }

    /**
     *
     * @param vc
     * @param refBasesStartingAtVCWithPad
     * @return
     */
    public static Pair<List<Integer>, byte[]> getNumTandemRepeatUnits(final VariantContext vc, final byte[] refBasesStartingAtVCWithPad) {
        Utils.nonNull(vc);
        Utils.nonNull(refBasesStartingAtVCWithPad);

        if ( ! vc.isIndel() ){ // only indels are tandem repeats
            return null;
        }
        final boolean VERBOSE = false;
        final String refBasesStartingAtVCWithoutPad = new String(refBasesStartingAtVCWithPad).substring(1);

        final Allele refAllele = vc.getReference();
        final byte[] refAlleleBases = Arrays.copyOfRange(refAllele.getBases(), 1, refAllele.length());

        byte[] repeatUnit = null;
        final List<Integer> lengths = new ArrayList<>();

        for ( final Allele allele : vc.getAlternateAlleles() ) {
            Pair<int[],byte[]> result = getNumTandemRepeatUnits(refAlleleBases, Arrays.copyOfRange(allele.getBases(), 1, allele.length()), refBasesStartingAtVCWithoutPad.getBytes());

            final int[] repetitionCount = result.getLeft();
            // repetition count = 0 means allele is not a tandem expansion of context
            if (repetitionCount[0] == 0 || repetitionCount[1] == 0)
                return null;

            if (lengths.isEmpty()) {
                lengths.add(repetitionCount[0]); // add ref allele length only once
            }
            lengths.add(repetitionCount[1]);  // add this alt allele's length

            repeatUnit = result.getRight();
            if (VERBOSE) {
                System.out.println("RefContext:"+refBasesStartingAtVCWithoutPad);
                System.out.println("Ref:"+refAllele.toString()+" Count:" + String.valueOf(repetitionCount[0]));
                System.out.println("Allele:"+allele.toString()+" Count:" + String.valueOf(repetitionCount[1]));
                System.out.println("RU:"+new String(repeatUnit));
            }
        }

        return new MutablePair<>(lengths,repeatUnit);
    }

    public static Pair<int[],byte[]> getNumTandemRepeatUnits(final byte[] refBases, final byte[] altBases, final byte[] remainingRefContext) {
         /* we can't exactly apply same logic as in basesAreRepeated() to compute tandem unit and number of repeated units.
           Consider case where ref =ATATAT and we have an insertion of ATAT. Natural description is (AT)3 -> (AT)2.
         */

        byte[] longB;
        // find first repeat unit based on either ref or alt, whichever is longer
        if (altBases.length > refBases.length)
            longB = altBases;
        else
            longB = refBases;

        // see if non-null allele (either ref or alt, whichever is longer) can be decomposed into several identical tandem units
        // for example, -*,CACA needs to first be decomposed into (CA)2
        final int repeatUnitLength = findRepeatedSubstring(longB);
        final byte[] repeatUnit = Arrays.copyOf(longB, repeatUnitLength);

        final int[] repetitionCount = new int[2];
        // look for repetitions forward on the ref bases (i.e. starting at beginning of ref bases)
        int repetitionsInRef = findNumberOfRepetitions(repeatUnit, refBases, true);
        repetitionCount[0] = findNumberOfRepetitions(repeatUnit, ArrayUtils.addAll(refBases, remainingRefContext), true)-repetitionsInRef;
        repetitionCount[1] = findNumberOfRepetitions(repeatUnit, ArrayUtils.addAll(altBases, remainingRefContext), true)-repetitionsInRef;

        return new MutablePair<>(repetitionCount, repeatUnit);

    }

    /**
     * Find out if a string can be represented as a tandem number of substrings.
     * For example ACTACT is a 2-tandem of ACT,
     * but ACTACA is not.
     *
     * @param bases                 String to be tested
     * @return                      Length of repeat unit, if string can be represented as tandem of substring (if it can't
     *                              be represented as one, it will be just the length of the input string)
     */
    public static int findRepeatedSubstring(byte[] bases) {

        int repLength;
        for (repLength=1; repLength <=bases.length; repLength++) {
            final byte[] candidateRepeatUnit = Arrays.copyOf(bases,repLength);
            boolean allBasesMatch = true;
            for (int start = repLength; start < bases.length; start += repLength ) {
                // check that remaining of string is exactly equal to repeat unit
                final byte[] basePiece = Arrays.copyOfRange(bases,start,start+candidateRepeatUnit.length);
                if (!Arrays.equals(candidateRepeatUnit, basePiece)) {
                    allBasesMatch = false;
                    break;
                }
            }
            if (allBasesMatch)
                return repLength;
        }

        return repLength;
    }

    /**
     * Finds number of repetitions a string consists of.
     * For example, for string ATAT and repeat unit AT, number of repetitions = 2
     * @param repeatUnit             Non-empty substring represented by byte array
     * @param testString             String to test (represented by byte array), may be empty
     * @param leadingRepeats         Look for leading (at the beginning of string) or trailing (at end of string) repetitions
     * For example:
     *    GATAT has 0 leading repeats of AT but 2 trailing repeats of AT
     *    ATATG has 1 leading repeat of A but 2 leading repeats of AT
     *    CCCCCCCC has 2 leading and 2 trailing repeats of CCC
     *
     * @return  Number of repetitions (0 if testString is not a concatenation of n repeatUnit's, including the case of empty testString)
     */
    public static int findNumberOfRepetitions(byte[] repeatUnit, byte[] testString, boolean leadingRepeats) {
        Utils.nonNull(repeatUnit, "repeatUnit");
        Utils.nonNull(testString, "testString");
        Utils.validateArg(repeatUnit.length != 0, "empty repeatUnit");
        if (testString.length == 0){
            return 0;
        }
        return findNumberOfRepetitions(repeatUnit, 0, repeatUnit.length, testString, 0, testString.length, leadingRepeats);
    }

    /**
     * Finds number of repetitions a string consists of.
     * Same as {@link #findNumberOfRepetitions} but operates on subarrays of a bigger array to save on copying.
     * For example, for string ATAT and repeat unit AT, number of repetitions = 2
     * @param repeatUnitFull             Non-empty substring represented by byte array
     * @param offsetInRepeatUnitFull     the offset in repeatUnitFull from which to read the repeat unit
     * @param repeatUnitLength           length of the repeat unit
     * @param testStringFull             string to test (represented by byte array), may be empty
     * @param offsetInTestStringFull     the offset in offsetInRepeatUnitFull from which to read the test string
     * @param testStringLength           length of the test string
     * @param leadingRepeats         Look for leading (at the beginning of string) or trailing (at end of string) repetitions
     * For example:
     *    GATAT has 0 leading repeats of AT but 2 trailing repeats of AT
     *    ATATG has 1 leading repeat of A but 2 leading repeats of AT
     *    CCCCCCCC has 2 leading and 2 trailing repeats of CCC
     * @return  Number of repetitions (0 if testString is not a concatenation of n repeatUnit's, including the case of empty testString)
     */
    public static int findNumberOfRepetitions(final byte[] repeatUnitFull, final int offsetInRepeatUnitFull, final int repeatUnitLength, final byte[] testStringFull, final int offsetInTestStringFull, final int testStringLength, final boolean leadingRepeats) {
        Utils.nonNull(repeatUnitFull, "repeatUnit");
        Utils.nonNull(testStringFull, "testString");
        Utils.validIndex(offsetInRepeatUnitFull, repeatUnitFull.length);
        Utils.validateArg(repeatUnitLength >= 0 && repeatUnitLength <= repeatUnitFull.length, "repeatUnitLength");
        if (testStringLength == 0){
            return 0;
        }
        Utils.validIndex(offsetInTestStringFull, testStringFull.length);
        Utils.validateArg(testStringLength >= 0 && testStringLength <= testStringFull.length, "testStringLength");
        final int lengthDifference = testStringLength - repeatUnitLength;

        if (leadingRepeats) {
            int numRepeats = 0;
            // look forward on the test string
            for (int start = 0; start <= lengthDifference; start += repeatUnitLength) {
                if(Utils.equalRange(testStringFull, start + offsetInTestStringFull, repeatUnitFull, offsetInRepeatUnitFull, repeatUnitLength)) {
                    numRepeats++;
                } else {
                    return numRepeats;
                }
            }
            return numRepeats;
        } else {
            // look backward. For example, if repeatUnit = AT and testString = GATAT, number of repeat units is still 2
            int numRepeats = 0;
            // look backward on the test string
            for (int start = lengthDifference; start >= 0; start -= repeatUnitLength) {
                if (Utils.equalRange(testStringFull, start + offsetInTestStringFull, repeatUnitFull, offsetInRepeatUnitFull, repeatUnitLength)) {
                    numRepeats++;
                } else {
                    return numRepeats;
                }
            }
            return numRepeats;
        }
    }

    /**
     * Helper function for isTandemRepeat that checks that allele matches somewhere on the reference
     */
    protected static boolean isRepeatAllele(final Allele ref, final Allele alt, final String refBasesStartingAtVCWithoutPad) {
        if ( ! Allele.oneIsPrefixOfOther(ref, alt) )
            return false; // we require one allele be a prefix of another

        if ( ref.length() > alt.length() ) { // we are a deletion
            return basesAreRepeated(ref.getBaseString(), alt.getBaseString(), refBasesStartingAtVCWithoutPad, 2);
        } else { // we are an insertion
            return basesAreRepeated(alt.getBaseString(), ref.getBaseString(), refBasesStartingAtVCWithoutPad, 1);
        }
    }

    protected static boolean basesAreRepeated(final String l, final String s, final String ref, final int minNumberOfMatches) {
        final String potentialRepeat = l.substring(s.length()); // skip s bases

        for ( int i = 0; i < minNumberOfMatches; i++) {
            final int start = i * potentialRepeat.length();
            final int end = (i+1) * potentialRepeat.length();
            if ( ref.length() < end )
                return false; // we ran out of bases to test
            final String refSub = ref.substring(start, end);
            if ( ! refSub.equals(potentialRepeat) )
                return false; // repeat didn't match, fail
        }

        return true; // we passed all tests, we matched
    }

    /**
     * Subset the samples in VC to reference only information with ref call alleles
     *
     * Preserves DP if present
     *
     * @param vc the variant context to subset down to
     * @param defaultPloidy defaultPloidy to use if a genotype doesn't have any alleles
     * @return a GenotypesContext
     */
    public static GenotypesContext subsetToRefOnly(final VariantContext vc, final int defaultPloidy) {
        Utils.nonNull(vc == null, "vc cannot be null");
        Utils.validateArg(defaultPloidy >= 1, () -> "defaultPloidy must be >= 1 but got " + defaultPloidy);

        final GenotypesContext oldGTs = vc.getGenotypes();
        if (oldGTs.isEmpty()) return oldGTs;
        final GenotypesContext newGTs = GenotypesContext.create(oldGTs.size());

        final Allele ref = vc.getReference();
        final List<Allele> diploidRefAlleles = Arrays.asList(ref, ref);

        // create the new genotypes
        for ( final Genotype g : vc.getGenotypes() ) {
            final int gPloidy = g.getPloidy() == 0 ? defaultPloidy : g.getPloidy();
            final List<Allele> refAlleles = gPloidy == 2 ? diploidRefAlleles : Collections.nCopies(gPloidy, ref);
            final GenotypeBuilder gb = new GenotypeBuilder(g.getSampleName(), refAlleles);
            if ( g.hasDP() ) gb.DP(g.getDP());
            if ( g.hasGQ() ) gb.GQ(g.getGQ());
            newGTs.add(gb.make());
        }

        return newGTs;
    }

    public static Genotype removePLsAndAD(final Genotype g) {
        return ( g.hasLikelihoods() || g.hasAD() ) ? new GenotypeBuilder(g).noPL().noAD().make() : g;
    }

    //TODO consider refactor variant-context merging code so that we share as much as possible between
    //TODO simpleMerge and referenceConfidenceMerge
    //TODO likely using a separate helper class or hierarchy.
    /**
     * Merges VariantContexts into a single hybrid.  Takes genotypes for common samples in priority order, if provided.
     * If uniquifySamples is true, the priority order is ignored and names are created by concatenating the VC name with
     * the sample name
     *
     * @param unsortedVCs               collection of unsorted VCs
     * @param priorityListOfVCs         priority list detailing the order in which we should grab the VCs
     * @param filteredRecordMergeType   merge type for filtered records
     * @param genotypeMergeOptions      merge option for genotypes
     * @param annotateOrigin            should we annotate the set it came from?
     * @param printMessages             should we print messages?
     * @param setKey                    the key name of the set
     * @param filteredAreUncalled       are filtered records uncalled?
     * @param mergeInfoWithMaxAC        should we merge in info from the VC with maximum allele count?
     * @return new VariantContext       representing the merge of unsortedVCs
     */
    public static VariantContext simpleMerge(final Collection<VariantContext> unsortedVCs,
                                             final List<String> priorityListOfVCs,
                                             final FilteredRecordMergeType filteredRecordMergeType,
                                             final GenotypeMergeType genotypeMergeOptions,
                                             final boolean annotateOrigin,
                                             final boolean printMessages,
                                             final String setKey,
                                             final boolean filteredAreUncalled,
                                             final boolean mergeInfoWithMaxAC ) {
        int originalNumOfVCs = priorityListOfVCs == null ? 0 : priorityListOfVCs.size();
        return simpleMerge(unsortedVCs, priorityListOfVCs, originalNumOfVCs, filteredRecordMergeType, genotypeMergeOptions, annotateOrigin, printMessages, setKey, filteredAreUncalled, mergeInfoWithMaxAC);
    }

    /**
     * Merges VariantContexts into a single hybrid.  Takes genotypes for common samples in priority order, if provided.
     * If uniquifySamples is true, the priority order is ignored and names are created by concatenating the VC name with
     * the sample name.
     * simpleMerge does not verify any more unique sample names EVEN if genotypeMergeOptions == GenotypeMergeType.REQUIRE_UNIQUE. One should use
     * SampleUtils.verifyUniqueSamplesNames to check that before using simpleMerge.
     *
     * For more information on this method see: http://www.thedistractionnetwork.com/programmer-problem/
     *
     * @param unsortedVCs               collection of unsorted VCs
     * @param priorityListOfVCs         priority list detailing the order in which we should grab the VCs
     * @param filteredRecordMergeType   merge type for filtered records
     * @param genotypeMergeOptions      merge option for genotypes
     * @param annotateOrigin            should we annotate the set it came from?
     * @param printMessages             should we print messages?
     * @param setKey                    the key name of the set
     * @param filteredAreUncalled       are filtered records uncalled?
     * @param mergeInfoWithMaxAC        should we merge in info from the VC with maximum allele count?
     * @return new VariantContext       representing the merge of unsortedVCs
     */
    public static VariantContext simpleMerge(final Collection<VariantContext> unsortedVCs,
                                             final List<String> priorityListOfVCs,
                                             final int originalNumOfVCs,
                                             final FilteredRecordMergeType filteredRecordMergeType,
                                             final GenotypeMergeType genotypeMergeOptions,
                                             final boolean annotateOrigin,
                                             final boolean printMessages,
                                             final String setKey,
                                             final boolean filteredAreUncalled,
                                             final boolean mergeInfoWithMaxAC ) {
        if ( unsortedVCs == null || unsortedVCs.isEmpty() )
            return null;

        if (priorityListOfVCs != null && originalNumOfVCs != priorityListOfVCs.size())
            throw new IllegalArgumentException("the number of the original VariantContexts must be the same as the number of VariantContexts in the priority list");

        if ( annotateOrigin && priorityListOfVCs == null && originalNumOfVCs == 0)
            throw new IllegalArgumentException("Cannot merge calls and annotate their origins without a complete priority list of VariantContexts or the number of original VariantContexts");

        final List<VariantContext> preFilteredVCs = sortVariantContextsByPriority(unsortedVCs, priorityListOfVCs, genotypeMergeOptions);
        // Make sure all variant contexts are padded with reference base in case of indels if necessary
        List<VariantContext> VCs = preFilteredVCs.stream()
                .filter(vc -> !filteredAreUncalled || vc.isNotFiltered())
                .collect(Collectors.toList());

        if ( VCs.isEmpty() ) // everything is filtered out and we're filteredAreUncalled
            return null;

        // establish the baseline info from the first VC
        final VariantContext first = VCs.get(0);
        final String name = first.getSource();
        final Allele refAllele = determineReferenceAllele(VCs);

        final LinkedHashSet<Allele> alleles = new LinkedHashSet<>();
        final Set<String> filters = new LinkedHashSet<>();
        final Map<String, Object> attributes = new LinkedHashMap<>();
        final Set<String> inconsistentAttributes = new LinkedHashSet<>();
        final Set<String> variantSources = new LinkedHashSet<>(); // contains the set of sources we found in our set of VCs that are variant
        final Set<String> rsIDs = new LinkedHashSet<>(1); // most of the time there's one id

        VariantContext longestVC = first;
        int depth = 0;
        int maxAC = -1;
        final Map<String, Object> attributesWithMaxAC = new LinkedHashMap<>();
        double log10PError = CommonInfo.NO_LOG10_PERROR;
        boolean anyVCHadFiltersApplied = false;
        VariantContext vcWithMaxAC = null;
        GenotypesContext genotypes = GenotypesContext.create();

        // counting the number of filtered and variant VCs
        int nFiltered = 0;

        boolean remapped = false;

        // cycle through and add info from the other VCs, making sure the loc/reference matches
        for ( final VariantContext vc : VCs ) {
            Utils.validate(longestVC.getStart() == vc.getStart(), () -> "BUG: attempting to merge VariantContexts with different start sites: first="+ first.toString() + " second=" + vc.toString());

            if ( VariantContextUtils.getSize(vc) > VariantContextUtils.getSize(longestVC) )
                longestVC = vc; // get the longest location

            nFiltered += vc.isFiltered() ? 1 : 0;
            if ( vc.isVariant() ) variantSources.add(vc.getSource());

            AlleleMapper alleleMapping = resolveIncompatibleAlleles(refAllele, vc, alleles);
            remapped = remapped || alleleMapping.needsRemapping();

            alleles.addAll(alleleMapping.values());

            mergeGenotypes(genotypes, vc, alleleMapping, genotypeMergeOptions == GenotypeMergeType.UNIQUIFY);

            // We always take the QUAL of the first VC with a non-MISSING qual for the combined value
            if ( log10PError == CommonInfo.NO_LOG10_PERROR )
                log10PError =  vc.getLog10PError();

            filters.addAll(vc.getFilters());
            anyVCHadFiltersApplied |= vc.filtersWereApplied();

            //
            // add attributes
            //
            // special case DP (add it up) and ID (just preserve it)
            //
            if (vc.hasAttribute(VCFConstants.DEPTH_KEY))
                depth += vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0);
            if ( vc.hasID() ) rsIDs.add(vc.getID());
            if (mergeInfoWithMaxAC && vc.hasAttribute(VCFConstants.ALLELE_COUNT_KEY)) {
                String rawAlleleCounts = vc.getAttributeAsString(VCFConstants.ALLELE_COUNT_KEY, null);
                // lets see if the string contains a "," separator
                if (rawAlleleCounts.contains(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)) {
                    final List<String> alleleCountArray = Arrays.asList(rawAlleleCounts.substring(1, rawAlleleCounts.length() - 1).split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR));
                    for (final String alleleCount : alleleCountArray) {
                        final int ac = Integer.valueOf(alleleCount.trim());
                        if (ac > maxAC) {
                            maxAC = ac;
                            vcWithMaxAC = vc;
                        }
                    }
                } else {
                    final int ac = Integer.valueOf(rawAlleleCounts);
                    if (ac > maxAC) {
                        maxAC = ac;
                        vcWithMaxAC = vc;
                    }
                }
            }

            for (final Map.Entry<String, Object> p : vc.getAttributes().entrySet()) {
                final String key = p.getKey();
                final Object value = p.getValue();
                // only output annotations that have the same value in every input VC
                // if we don't like the key already, don't go anywhere
                if ( ! inconsistentAttributes.contains(key) ) {
                    final boolean alreadyFound = attributes.containsKey(key);
                    final Object boundValue = attributes.get(key);
                    final boolean boundIsMissingValue = alreadyFound && boundValue.equals(VCFConstants.MISSING_VALUE_v4);

                    if ( alreadyFound && ! boundValue.equals(value) && ! boundIsMissingValue ) {
                        // we found the value but we're inconsistent, put it in the exclude list
                        inconsistentAttributes.add(key);
                        attributes.remove(key);
                    } else if ( ! alreadyFound || boundIsMissingValue )  { // no value
                        attributes.put(key, value);
                    }
                }
            }
        }

        // if we have more alternate alleles in the merged VC than in one or more of the
        // original VCs, we need to strip out the GL/PLs (because they are no longer accurate), as well as allele-dependent attributes like AC,AF, and AD
        for ( final VariantContext vc : VCs ) {
            if (vc.getAlleles().size() == 1)
                continue;
            if ( hasPLIncompatibleAlleles(alleles, vc.getAlleles())) {
                if ( ! genotypes.isEmpty() ) {
                    logger.debug(String.format("Stripping PLs at %s:%d-%d due to incompatible alleles merged=%s vs. single=%s",
                            vc.getContig(), vc.getStart(), vc.getEnd(), alleles, vc.getAlleles()));
                }
                genotypes = stripPLsAndAD(genotypes);
                // this will remove stale AC,AF attributed from vc
                VariantContextUtils.calculateChromosomeCounts(vc, attributes, true);
                break;
            }
        }

        // take the VC with the maxAC and pull the attributes into a modifiable map
        if ( mergeInfoWithMaxAC && vcWithMaxAC != null ) {
            attributesWithMaxAC.putAll(vcWithMaxAC.getAttributes());
        }

        // if at least one record was unfiltered and we want a union, clear all of the filters
        if ( (filteredRecordMergeType == FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED && nFiltered != VCs.size()) || filteredRecordMergeType == FilteredRecordMergeType.KEEP_UNCONDITIONAL )
            filters.clear();


        if ( annotateOrigin ) { // we care about where the call came from
            String setValue;
            if ( nFiltered == 0 && variantSources.size() == originalNumOfVCs ) // nothing was unfiltered
                setValue = MERGE_INTERSECTION;
            else if ( nFiltered == VCs.size() )     // everything was filtered out
                setValue = MERGE_FILTER_IN_ALL;
            else if ( variantSources.isEmpty() )    // everyone was reference
                setValue = MERGE_REF_IN_ALL;
            else {
                final LinkedHashSet<String> s = new LinkedHashSet<>();
                for ( final VariantContext vc : VCs )
                    if ( vc.isVariant() )
                        s.add( vc.isFiltered() ? MERGE_FILTER_PREFIX + vc.getSource() : vc.getSource() );
                setValue = Utils.join("-", s);
            }

            if ( setKey != null ) {
                attributes.put(setKey, setValue);
                if( mergeInfoWithMaxAC && vcWithMaxAC != null ) {
                    attributesWithMaxAC.put(setKey, setValue);
                }
            }
        }

        if ( depth > 0 )
            attributes.put(VCFConstants.DEPTH_KEY, String.valueOf(depth));

        final String ID = rsIDs.isEmpty() ? VCFConstants.EMPTY_ID_FIELD : Utils.join(",", rsIDs);

        final VariantContextBuilder builder = new VariantContextBuilder().source(name).id(ID);
        builder.loc(longestVC.getContig(), longestVC.getStart(), longestVC.getEnd());
        builder.alleles(alleles);
        builder.genotypes(genotypes);
        builder.log10PError(log10PError);
        if ( anyVCHadFiltersApplied ) {
            builder.filters(filters.isEmpty() ? filters : new TreeSet<>(filters));
        }
        builder.attributes(new TreeMap<>(mergeInfoWithMaxAC ? attributesWithMaxAC : attributes));

        // Trim the padded bases of all alleles if necessary
        final VariantContext merged = builder.make();
        if ( printMessages && remapped ) System.out.printf("Remapped => %s%n", merged);
        return merged;
    }


    //TODO as part of a larger refactoring effort remapAlleles can be merged with createAlleleMapping.

    public static GenotypesContext stripPLsAndAD(final GenotypesContext genotypes) {
        final GenotypesContext newGs = GenotypesContext.create(genotypes.size());
        for ( final Genotype g : genotypes ) {
            newGs.add(removePLsAndAD(g));
        }
        return newGs;
    }

    private static Allele determineReferenceAllele(final List<VariantContext> VCs) {
        return determineReferenceAllele(VCs, null);
    }

    public static boolean contextMatchesLoc(final VariantContext vc, final Locatable loc) {
        return loc == null || loc.getStart() == vc.getStart();
    }

    public static AlleleMapper resolveIncompatibleAlleles(final Allele refAllele, final VariantContext vc, final LinkedHashSet<Allele> allAlleles) {
        if ( refAllele.equals(vc.getReference()) )
            return new AlleleMapper(vc);
        else {
            final Map<Allele, Allele> map = createAlleleMapping(refAllele, vc, allAlleles);
            map.put(vc.getReference(), refAllele);
            return new AlleleMapper(map);
        }
    }

    //TODO as part of a larger refactoring effort {@link #createAlleleMapping} can be merged with {@link ReferenceConfidenceVariantContextMerger#remapAlleles}.
    /**
     * Create an allele mapping for the given context where its reference allele must (potentially) be extended to the given allele
     *
     * The refAllele is the longest reference allele seen at this start site.
     * So imagine it is:
     * refAllele: ACGTGA
     * myRef:     ACGT
     * myAlt:     A
     *
     * We need to remap all of the alleles in vc to include the extra GA so that
     * myRef => refAllele and myAlt => AGA
     *
     * @param refAllele          the new (extended) reference allele
     * @param oneVC              the Variant Context to extend
     * @param currentAlleles     the list of alleles already created
     * @return a non-null mapping of original alleles to new (extended) ones
     */
    public static Map<Allele, Allele> createAlleleMapping(final Allele refAllele,
                                                           final VariantContext oneVC,
                                                           final Collection<Allele> currentAlleles) {
        final Allele myRef = oneVC.getReference();
        Utils.validate(refAllele.length() > myRef.length(), () -> "BUG: myRef="+myRef+" is longer than refAllele="+refAllele);
        final byte[] extraBases = Arrays.copyOfRange(refAllele.getBases(), myRef.length(), refAllele.length());

        final Map<Allele, Allele> map = new LinkedHashMap<>();
        for ( final Allele a : oneVC.getAlternateAlleles() ) {
            if ( isNonSymbolicExtendableAllele(a) ) {
                Allele extended = Allele.extend(a, extraBases);
                for ( final Allele b : currentAlleles )
                    if ( extended.equals(b) )
                        extended = b;
                map.put(a, extended);
            } else if (a.equals(Allele.SPAN_DEL)) {
                map.put(a, a);
            }
        }

        return map;
    }

    private static boolean isNonSymbolicExtendableAllele(final Allele allele) {
        return ! (allele.isReference() || allele.isSymbolic() || allele.equals(Allele.SPAN_DEL));
    }

    public static List<VariantContext> sortVariantContextsByPriority(Collection<VariantContext> unsortedVCs, List<String> priorityListOfVCs, GenotypeMergeType mergeOption ) {
        if ( mergeOption == GenotypeMergeType.PRIORITIZE && priorityListOfVCs == null )
            throw new IllegalArgumentException("Cannot merge calls by priority with a null priority list");

        if ( priorityListOfVCs == null || mergeOption == GenotypeMergeType.UNSORTED )
            return new ArrayList<>(unsortedVCs);
        else {
            ArrayList<VariantContext> sorted = new ArrayList<>(unsortedVCs);
            Collections.sort(sorted, new CompareByPriority(priorityListOfVCs));
            return sorted;
        }
    }

    private static void mergeGenotypes(GenotypesContext mergedGenotypes, VariantContext oneVC, AlleleMapper alleleMapping, boolean uniquifySamples) {
        //TODO: should we add a check for cases when the genotypeMergeOption is REQUIRE_UNIQUE
        for ( final Genotype g : oneVC.getGenotypes() ) {
            final String name = mergedSampleName(oneVC.getSource(), g.getSampleName(), uniquifySamples);
            if ( ! mergedGenotypes.containsSample(name) ) {
                // only add if the name is new
                Genotype newG = g;

                if ( uniquifySamples || alleleMapping.needsRemapping() ) {
                    final List<Allele> alleles = alleleMapping.needsRemapping() ? alleleMapping.remap(g.getAlleles()) : g.getAlleles();
                    newG = new GenotypeBuilder(g).name(name).alleles(alleles).make();
                }

                mergedGenotypes.add(newG);
            }
        }
    }

    /**
     * Cached NO_CALL immutable lists where the position ith contains the list with i elements.
     */
    @SuppressWarnings({"rawtypes", "unchecked"})
    private static List<Allele>[] NOCALL_LISTS = new List[] {
            Collections.emptyList(),
            Collections.singletonList(Allele.NO_CALL),
            Collections.nCopies(2,Allele.NO_CALL)
    };

    /**
     * Code to ensure that {@link #NOCALL_LISTS} has enough entries beyond the requested ploidy
     * @param capacity the requested ploidy.
     */
    private static void ensureNoCallListsCapacity(final int capacity) {
        final int currentCapacity = NOCALL_LISTS.length - 1;
        if (currentCapacity >= capacity)
            return;
        NOCALL_LISTS = Arrays.copyOf(NOCALL_LISTS,Math.max(capacity,currentCapacity << 1) + 1);
        for (int i = currentCapacity + 1; i < NOCALL_LISTS.length; i++)
            NOCALL_LISTS[i] = Collections.nCopies(i,Allele.NO_CALL);
    }

    /**
     * Returns a {@link Allele#NO_CALL NO_CALL} allele list provided the ploidy.
     *
     * @param ploidy the required ploidy.
     *
     * @return never {@code null}, but an empty list if {@code ploidy} is equal or less than 0. The returned list
     *   might or might not be mutable.
     */
    public static List<Allele> noCallAlleles(final int ploidy) {
        if (NOCALL_LISTS.length <= ploidy)
            ensureNoCallListsCapacity(ploidy);
        return NOCALL_LISTS[ploidy];
    }


    public static String mergedSampleName(String trackName, String sampleName, boolean uniquify ) {
        return uniquify ? sampleName + "." + trackName : sampleName;
    }

    /**
     * Trim the alleles in inputVC from the reverse direction
     *
     * @param inputVC a non-null input VC whose alleles might need a haircut
     * @return a non-null VariantContext (may be == to inputVC) with alleles trimmed up
     */
    public static VariantContext reverseTrimAlleles( final VariantContext inputVC ) {
        return trimAlleles(inputVC, false, true);
    }

    /**
     * Trim the alleles in inputVC forward and reverse, as requested
     *
     * @param inputVC a non-null input VC whose alleles might need a haircut
     * @param trimForward should we trim up the alleles from the forward direction?
     * @param trimReverse should we trim up the alleles from the reverse direction?
     * @return a non-null VariantContext (may be == to inputVC) with trimmed up alleles
     */
    public static VariantContext trimAlleles(final VariantContext inputVC, final boolean trimForward, final boolean trimReverse) {
        Utils.nonNull(inputVC);

        if ( inputVC.getNAlleles() <= 1 || inputVC.isSNP() )
            return inputVC;

        // see whether we need to trim common reference base from all alleles
        final int revTrim = trimReverse ? computeReverseClipping(inputVC.getAlleles(), inputVC.getReference().getDisplayString().getBytes()) : 0;
        final VariantContext revTrimVC = trimAlleles(inputVC, -1, revTrim);
        final int fwdTrim = trimForward ? computeForwardClipping(revTrimVC.getAlleles()) : -1;
        return trimAlleles(revTrimVC, fwdTrim, 0);
    }

    /**
     * Trim up alleles in inputVC, cutting out all bases up to fwdTrimEnd inclusive and
     * the last revTrim bases from the end
     *
     * @param inputVC a non-null input VC
     * @param fwdTrimEnd bases up to this index (can be -1) will be removed from the start of all alleles
     * @param revTrim the last revTrim bases of each allele will be clipped off as well
     * @return a non-null VariantContext (may be == to inputVC) with trimmed up alleles
     */
    protected static VariantContext trimAlleles(final VariantContext inputVC,
                                                final int fwdTrimEnd,
                                                final int revTrim) {
        if( fwdTrimEnd == -1 && revTrim == 0 ) // nothing to do, so just return inputVC unmodified
            return inputVC;

        final List<Allele> alleles = new LinkedList<>();
        final Map<Allele, Allele> originalToTrimmedAlleleMap = new LinkedHashMap<>();

        for (final Allele a : inputVC.getAlleles()) {
            if (a.isSymbolic()) {
                alleles.add(a);
                originalToTrimmedAlleleMap.put(a, a);
            } else {
                // get bases for current allele and create a new one with trimmed bases
                final byte[] newBases = Arrays.copyOfRange(a.getBases(), fwdTrimEnd+1, a.length()-revTrim);
                final Allele trimmedAllele = Allele.create(newBases, a.isReference());
                alleles.add(trimmedAllele);
                originalToTrimmedAlleleMap.put(a, trimmedAllele);
            }
        }

        // now we can recreate new genotypes with trimmed alleles
        final AlleleMapper alleleMapper = new AlleleMapper(originalToTrimmedAlleleMap);
        final GenotypesContext genotypes = updateGenotypesWithMappedAlleles(inputVC.getGenotypes(), alleleMapper);

        final int start = inputVC.getStart() + (fwdTrimEnd + 1);
        final VariantContextBuilder builder = new VariantContextBuilder(inputVC);
        builder.start(start);
        builder.stop(start + alleles.get(0).length() - 1);
        builder.alleles(alleles);
        builder.genotypes(genotypes);
        return builder.make();
    }

    protected static GenotypesContext updateGenotypesWithMappedAlleles(final GenotypesContext originalGenotypes, final AlleleMapper alleleMapper) {
        final GenotypesContext updatedGenotypes = GenotypesContext.create(originalGenotypes.size());

        for ( final Genotype genotype : originalGenotypes ) {
            final List<Allele> updatedAlleles = alleleMapper.remap(genotype.getAlleles());
            updatedGenotypes.add(new GenotypeBuilder(genotype).alleles(updatedAlleles).make());
        }

        return updatedGenotypes;
    }

    public static int computeReverseClipping(final List<Allele> unclippedAlleles, final byte[] ref) {
        int clipping = 0;
        boolean stillClipping = true;

        while ( stillClipping ) {
            for ( final Allele a : unclippedAlleles ) {
                if ( a.isSymbolic() )
                    continue;

                // we need to ensure that we don't reverse clip out all of the bases from an allele because we then will have the wrong
                // position set for the VariantContext (although it's okay to forward clip it all out, because the position will be fine).
                if ( a.length() - clipping == 0 )
                    return clipping - 1;

                if ( a.length() - clipping <= 0 || a.length() == 0 ) {
                    stillClipping = false;
                }
                else if ( ref.length == clipping ) {
                    return -1;
                }
                else if ( a.getBases()[a.length()-clipping-1] != ref[ref.length-clipping-1] ) {
                    stillClipping = false;
                }
            }
            if ( stillClipping )
                clipping++;
        }

        return clipping;
    }

    /**
     * Clip out any unnecessary bases off the front of the alleles
     *
     * The VCF spec represents alleles as block substitutions, replacing AC with A for a
     * 1 bp deletion of the C.  However, it's possible that we'd end up with alleles that
     * contain extra bases on the left, such as GAC/GA to represent the same 1 bp deletion.
     * This routine finds an offset among all alleles that can be safely trimmed
     * off the left of each allele and still represent the same block substitution.
     *
     * A/C => A/C
     * AC/A => AC/A
     * ACC/AC => CC/C
     * AGT/CAT => AGT/CAT
     * <DEL>/C => <DEL>/C
     *
     * @param unclippedAlleles a non-null list of alleles that we want to clip
     * @return the offset into the alleles where we can safely clip, inclusive, or
     *   -1 if no clipping is tolerated.  So, if the result is 0, then we can remove
     *   the first base of every allele.  If the result is 1, we can remove the
     *   second base.
     */
    public static int computeForwardClipping(final List<Allele> unclippedAlleles) {
        // cannot clip unless there's at least 1 alt allele
        if ( unclippedAlleles.size() <= 1 )
            return -1;

        // we cannot forward clip any set of alleles containing a symbolic allele
        int minAlleleLength = Integer.MAX_VALUE;
        for ( final Allele a : unclippedAlleles ) {
            if ( a.isSymbolic() )
                return -1;
            minAlleleLength = Math.min(minAlleleLength, a.length());
        }

        final byte[] firstAlleleBases = unclippedAlleles.get(0).getBases();
        int indexOflastSharedBase = -1;

        // the -1 to the stop is that we can never clip off the right most base
        for ( int i = 0; i < minAlleleLength - 1; i++) {
            final byte base = firstAlleleBases[i];

            for ( final Allele allele : unclippedAlleles ) {
                if ( allele.getBases()[i] != base )
                    return indexOflastSharedBase;
            }

            indexOflastSharedBase = i;
        }

        return indexOflastSharedBase;
    }

    protected static class AlleleMapper {
        private VariantContext vc = null;
        private Map<Allele, Allele> map = null;
        public AlleleMapper(VariantContext vc)          { this.vc = vc; }
        public AlleleMapper(Map<Allele, Allele> map)    { this.map = map; }
        public boolean needsRemapping()                 { return this.map != null; }
        public Collection<Allele> values()              { return map != null ? map.values() : vc.getAlleles(); }
        public Allele remap(Allele a)                   { return map != null && map.containsKey(a) ? map.get(a) : a; }

        public List<Allele> remap(List<Allele> as) {
            List<Allele> newAs = as.stream()
                    .map(this::remap)
                    .collect(Collectors.toList());
            return newAs;
        }

    }

    private static class CompareByPriority implements Comparator<VariantContext>, Serializable {
        private static final long serialVersionUID = 0L;

        List<String> priorityListOfVCs;
        public CompareByPriority(List<String> priorityListOfVCs) {
            this.priorityListOfVCs = priorityListOfVCs;
        }

        private int getIndex(VariantContext vc) {
            return Utils.validIndex(priorityListOfVCs.indexOf(vc.getSource()), priorityListOfVCs.size());
        }

        @Override
        public int compare(VariantContext vc1, VariantContext vc2) {
            return Integer.valueOf(getIndex(vc1)).compareTo(getIndex(vc2));
        }
    }

    /**
     * For testing purposes only.  Create a site-only VariantContext at contig:start containing alleles
     *
     * @param name the name of the VC
     * @param contig the contig for the VC
     * @param start the start of the VC
     * @param alleleStrings a non-null, non-empty list of strings for the alleles.  The first will be the ref allele, and others the
     *                      alt.  Will compute the stop of the VC from the length of the reference allele
     * @return a non-null VariantContext
     */
    public static VariantContext makeFromAlleles(final String name, final String contig, final int start, final List<String> alleleStrings) {
        Utils.nonNull(alleleStrings, "alleleStrings is null");
        Utils.nonEmpty(alleleStrings, "alleleStrings is empty");

        final List<Allele> alleles = new LinkedList<>();
        final int length = alleleStrings.get(0).length();

        boolean first = true;
        for ( final String alleleString : alleleStrings ) {
            alleles.add(Allele.create(alleleString, first));
            first = false;
        }
        return new VariantContextBuilder(name, contig, start, start+length-1, alleles).make();
    }

    /**
     * Split variant context into its biallelic components if there are more than 2 alleles
     * <p>
     * For VC has A/B/C alleles, returns A/B and A/C contexts.
     * Alleles are right trimmed to satisfy VCF conventions
     * <p>
     * If vc is biallelic or non-variant it is just returned
     * <p>
     * Chromosome counts are updated (but they are by definition 0)
     *
     * @param vc                       a potentially multi-allelic variant context
     * @param trimLeft                 if true, we will also left trim alleles, potentially moving the resulting vcs forward on the genome
     * @param genotypeAssignmentMethod assignment strategy for the (subsetted) PLs
     * @param keepOriginalChrCounts    keep the orignal chromosome counts before subsetting
     * @return a list of bi-allelic (or monomorphic) variant context
     */
    public static List<VariantContext> splitVariantContextToBiallelics(final VariantContext vc, final boolean trimLeft, final GenotypeAssignmentMethod genotypeAssignmentMethod,
                                                                       final boolean keepOriginalChrCounts) {
        Utils.nonNull(vc);

        if (!vc.isVariant() || vc.isBiallelic())
            // non variant or biallelics already satisfy the contract
            return Collections.singletonList(vc);
        else {
            final List<VariantContext> biallelics = new LinkedList<>();

            // if any of the genotypes are het-not-ref (i.e. 1/2), set all of them to no-call
            final GenotypeAssignmentMethod genotypeAssignmentMethodUsed = hasHetNonRef(vc.getGenotypes()) ? GenotypeAssignmentMethod.SET_TO_NO_CALL_NO_ANNOTATIONS : genotypeAssignmentMethod;

            for (final Allele alt : vc.getAlternateAlleles()) {
                final VariantContextBuilder builder = new VariantContextBuilder(vc);

                // make biallelic alleles
                final List<Allele> alleles = Arrays.asList(vc.getReference(), alt);
                builder.alleles(alleles);

                // since the VC has been subset, remove the invalid attributes
                for (final String key : vc.getAttributes().keySet()) {
                    if (!(key.equals(VCFConstants.ALLELE_COUNT_KEY) || key.equals(VCFConstants.ALLELE_FREQUENCY_KEY) || key.equals(VCFConstants.ALLELE_NUMBER_KEY)) ||
                            genotypeAssignmentMethodUsed == GenotypeAssignmentMethod.SET_TO_NO_CALL_NO_ANNOTATIONS) {
                        builder.rmAttribute(key);
                    }
                }

                // subset INFO field annotations if available if genotype is called
                if (genotypeAssignmentMethodUsed != GenotypeAssignmentMethod.SET_TO_NO_CALL_NO_ANNOTATIONS &&
                        genotypeAssignmentMethodUsed != GenotypeAssignmentMethod.SET_TO_NO_CALL)
                    AlleleSubsettingUtils.addInfoFieldAnnotations(vc, builder, keepOriginalChrCounts);

                builder.genotypes(AlleleSubsettingUtils.subsetAlleles(vc.getGenotypes(),2,vc.getAlleles(), alleles, genotypeAssignmentMethodUsed,vc.getAttributeAsInt("DP",0)));
                final VariantContext trimmed = trimAlleles(builder.make(), trimLeft, true);
                biallelics.add(trimmed);
            }

            return biallelics;
        }
    }

    /**
     * Check if any of the genotypes is heterozygous, non-reference (i.e. 1/2)
     *
     * @param genotypesContext genotype information
     * @return true if any of the genotypes are heterozygous, non-reference, false otherwise
     */
    private static boolean hasHetNonRef(final GenotypesContext genotypesContext) {
        for (final Genotype gt : genotypesContext) {
            if (gt.isHetNonRef())
                return true;
        }
        return false;
    }

    /**
     * Splits the alleles for the provided variant context into its primitive parts.
     * Requires that the input VC be bi-allelic, so calling methods should first call splitVariantContextToBiallelics() if needed.
     * Currently works only for MNPs.
     *
     * @param vc  the non-null VC to split
     * @return a non-empty list of VCs split into primitive parts or the original VC otherwise
     */
    public static List<VariantContext> splitIntoPrimitiveAlleles(final VariantContext vc) {
        Utils.nonNull(vc);
        if ( !vc.isBiallelic() ) {
            throw new IllegalArgumentException("Trying to break a multi-allelic Variant Context into primitive parts");
        }

        // currently only works for MNPs
        if ( !vc.isMNP() )
            return Arrays.asList(vc);

        final byte[] ref = vc.getReference().getBases();
        final byte[] alt = vc.getAlternateAllele(0).getBases();

        Utils.validate(ref.length == alt.length, "ref and alt alleles for MNP have different lengths");

        final List<VariantContext> result = new ArrayList<>(ref.length);

        for ( int i = 0; i < ref.length; i++ ) {

            // if the ref and alt bases are different at a given position, create a new SNP record (otherwise do nothing)
            if ( ref[i] != alt[i] ) {

                // create the ref and alt SNP alleles
                final Allele newRefAllele = Allele.create(ref[i], true);
                final Allele newAltAllele = Allele.create(alt[i], false);

                // create a new VariantContext with the new SNP alleles
                final VariantContextBuilder newVC = new VariantContextBuilder(vc).start(vc.getStart() + i).stop(vc.getStart() + i).alleles(Arrays.asList(newRefAllele, newAltAllele));

                // create new genotypes with updated alleles
                final Map<Allele, Allele> alleleMap = new LinkedHashMap<>();
                alleleMap.put(vc.getReference(), newRefAllele);
                alleleMap.put(vc.getAlternateAllele(0), newAltAllele);
                final GenotypesContext newGenotypes = updateGenotypesWithMappedAlleles(vc.getGenotypes(), new AlleleMapper(alleleMap));

                result.add(newVC.genotypes(newGenotypes).make());
            }
        }

        if ( result.isEmpty() )
            result.add(vc);

        return result;
    }

    /**
     * Are vc1 and 2 equal including their position and alleles?
     * @param vc1 non-null VariantContext
     * @param vc2 non-null VariantContext
     * @return true if vc1 and vc2 are equal, false otherwise
     */
    public static boolean equalSites(final VariantContext vc1, final VariantContext vc2) {
        Utils.nonNull(vc1, "vc1 is null");
        Utils.nonNull(vc2, "vc2 is null");

        if ( vc1.getStart() != vc2.getStart() ) return false;
        if ( vc1.getEnd() != vc2.getEnd() ) return false;
        if ( !vc1.getContig().equals(vc2.getContig())) return false;
        return vc1.getAlleles().equals(vc2.getAlleles());
    }

    /**
     * Returns the absolute 0-based index of an allele.
     *
     * <p/>
     * If the allele is equal to the reference, the result is 0, if it equal to the first alternative the result is 1
     * and so forth.
     * <p/>
     * Therefore if you want the 0-based index within the alternative alleles you need to do the following:
     *
     * <p/>
     * You can indicate whether the Java object reference comparator {@code ==} can be safelly used by setting {@code useEquals} to {@code false}.
     *
     * @param vc the target variant context.
     * @param allele the target allele.
     * @param ignoreRefState whether the reference states of the allele is important at all. Has no effect if {@code useEquals} is {@code false}.
     * @param considerRefAllele whether the reference allele should be considered. You should set it to {@code false} if you are only interested in alternative alleles.
     * @param useEquals whether equal method should be used in the search: {@link Allele#equals(Allele,boolean)}.
     *
     * @throws IllegalArgumentException if {@code allele} is {@code null}.
     * @return {@code -1} if there is no such allele that satify those criteria, a value between 0 and {@link VariantContext#getNAlleles()} {@code -1} otherwise.
     */
    public static int indexOfAllele(final VariantContext vc, final Allele allele, final boolean ignoreRefState, final boolean considerRefAllele, final boolean useEquals) {
        Utils.nonNull(allele);
        return useEquals ? indexOfEqualAllele(vc,allele,ignoreRefState,considerRefAllele) : indexOfSameAllele(vc,allele,considerRefAllele);
    }

    /**
     * Returns the relative 0-based index of an alternative allele.
     * <p/>
     * The the query allele is the same as the first alternative allele, the result is 0,
     * if it is equal to the second 1 and so forth.
     *
     *
     * <p/>
     * Notice that the ref-status of the query {@code allele} is ignored.
     *
     * @param vc the target variant context.
     * @param allele the query allele.
     * @param useEquals  whether equal method should be used in the search: {@link Allele#equals(Allele,boolean)}.
     *
     * @throws IllegalArgumentException if {@code allele} is {@code null}.
     *
     * @return {@code -1} if there is no such allele that satisfy those criteria, a value between 0 and the number
     *  of alternative alleles - 1.
     */
    public static int indexOfAltAllele(final VariantContext vc, final Allele allele, final boolean useEquals) {
        final int absoluteIndex = indexOfAllele(vc,allele,true,false,useEquals);
        return absoluteIndex == -1 ? -1 : absoluteIndex - 1;
    }

    // Implements index search using equals.
    private static int indexOfEqualAllele(final VariantContext vc, final Allele allele, final boolean ignoreRefState,
                                          final boolean considerRefAllele) {
        int i = 0;
        for (final Allele a : vc.getAlleles())
            if (a.equals(allele,ignoreRefState))
                return i == 0 ? (considerRefAllele ? 0 : -1) : i;
            else
                i++;
        return -1;
    }

    // Implements index search using ==.
    private static int indexOfSameAllele(final VariantContext vc, final Allele allele, final boolean considerRefAllele) {
        int i = 0;

        for (final Allele a : vc.getAlleles())
            if (a == allele)
                return i == 0 ? (considerRefAllele ? 0 : -1) : i;
            else
                i++;

        return -1;
    }

    /**
     * Add chromosome counts (AC, AN and AF) to the VCF header lines
     *
     * @param headerLines the VCF header lines
     */
    public static void addChromosomeCountsToHeader(final Set<VCFHeaderLine> headerLines) {
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY));
    }

    /**
     * Set the builder's filtered genotypes to no-call and update AC, AN and AF
     *
     * @param builder builder for variant context
     * @param vc the VariantContext record to set filtered genotypes to no-call
     * @param setFilteredGenotypesToNocall flag to set filtered genotype to NO CALL
     * @param filters the filters for each genotype
     */
    public static void setFilteredGenotypeToNocall(final VariantContextBuilder builder, final VariantContext vc,
                                                   final boolean setFilteredGenotypesToNocall, BiFunction<VariantContext, Genotype, List<String>> filters) {
        Utils.nonNull(vc);
        Utils.nonNull(builder);
        Utils.nonNull(filters);

        final GenotypesContext genotypes = GenotypesContext.create(vc.getGenotypes().size());

        // update genotypes if filtered genotypes are set to no-call
        boolean haveFilteredNoCallAlleles = false;
        for (final Genotype g : vc.getGenotypes()) {
            if (g.isCalled()) {
                List<String> filterNames = filters.apply(vc,g);
                if (!filterNames.isEmpty() && setFilteredGenotypesToNocall) {
                    haveFilteredNoCallAlleles = true;
                    genotypes.add(new GenotypeBuilder(g).filters(filterNames).alleles(GATKVariantContextUtils.noCallAlleles((g.getPloidy()))).make());
                } else {
                    genotypes.add(new GenotypeBuilder(g).filters(filterNames).make());
                }
            } else {
                genotypes.add(g);
            }
        }

        // if filtered genotypes are set to no-call, output recomputed AC, AN, AF
        if ( haveFilteredNoCallAlleles ) {
            final Map<String, Object> attributes = new LinkedHashMap<>(vc.getAttributes()); // need to make mutable
            VariantContextUtils.calculateChromosomeCounts(builder.genotypes(genotypes).make(), attributes, true, vc.getSampleNames());
            builder.attributes(attributes);
        }

        // update genotypes
        builder.genotypes(genotypes);
    }

    /**
     * Calculate the Genotype Quality (GQ) by subtracting the smallest Phred-scaled likelihoods (PL) from the second smallest.
     *
     * @param plValues  array of PL values
     * @return          the genotype quality corresponding to the PL values
     */
    public static int calculateGQFromPLs(final int[] plValues) {
        Utils.nonNull(plValues);
        Utils.validateArg(plValues.length >= 2, () -> "Array of PL values must contain at least two elements.");

        int first = plValues[0];
        int second = plValues[1];
        if (first > second) {
            second = first;
            first = plValues[1];
        }
        for (int i = 2; i < plValues.length; i++) {
            final int candidate = plValues[i];
            if (candidate >= second) {
                continue;
            }
            if (candidate <= first) {
                second = first;
                first = candidate;
            } else
                second = candidate;
        }
        return second - first;
    }

    /**
     *  Create a string with existing filters plus the one to append
     * @param vc VariantContext to base initial filter values.  Not {@code null}
     * @param filterToAppend the filter value to append to the strings Not {@code null}
     * @return String of filter values.  Sorted alphabetically.
     */
    public static List<String> createFilterListWithAppend(VariantContext vc, String filterToAppend) {
        Utils.nonNull(vc);
        Utils.nonNull(filterToAppend);

        final List<String> filtersAsList = new ArrayList<>(vc.getFilters());
        filtersAsList.add(filterToAppend);
        filtersAsList.sort(String::compareToIgnoreCase);
        return filtersAsList;
    }

    /**
     * Determines whether the given reference and alternate alleles constitute a frameshift mutation.
     * @param reference The reference {@link Allele}.  Must not be {@code null}.
     * @param alternate The alternate / variant {@link Allele}.  Must not be {@code null}.
     * @return {@code true} if replacing the reference with the alternate results in a frameshift.  {@code false} otherwise.
     */
    public static boolean isFrameshift(final Allele reference, final Allele alternate) {

        Utils.nonNull(reference);
        Utils.nonNull(alternate);

        // We know it's a frameshift if we have a replacement that is not of a
        // length evenly divisible by 3 because that's how many bases are read at once:
        return ((Math.abs( reference.length() - alternate.length() ) % 3) != 0);
    }

    /**
     * Determines whether the given reference and alternate alleles constitute a frameshift mutation.
     * @param reference The {@link String} containing the bases of the reference allele.  Must not be {@code null}.
     * @param alternate The {@link String} containing the bases of the alternate allele.  Must not be {@code null}.
     * @return {@code true} if replacing the reference with the alternate results in a frameshift.  {@code false} otherwise.
     */
    public static boolean isFrameshift(final String reference, final String alternate) {

        Utils.nonNull(reference);
        Utils.nonNull(alternate);

        final String refComparable = reference.replaceAll("\\*", "");
        final String altComperable = alternate.replaceAll("\\*", "");

        // We know it's a frameshift if we have a replacement that is not of a
        // length evenly divisible by 3 because that's how many bases are read at once:
        return ((Math.abs( refComparable.length() - altComperable.length() ) % 3) != 0);
    }

    /**
     * Determines whether the given reference and alternate alleles constitute a frameshift mutation.
     * This does not take into account the cases where either or both of the alleles overlap splice sites.
     * @param startPos Genomic start position (1-based, inclusive) of the variant.  Must be > 0.
     * @param refEnd Genomic end position (1-based, inclusive) of the reference allele.  Must be > 0.
     * @param altEnd Genomic end position (1-based, inclusive) of the alternate allele.  Must be > 0.
     * @return {@code true} if replacing the reference with the alternate results in a frameshift.  {@code false} otherwise.
     */
    public static boolean isFrameshift(final int startPos, final int refEnd, final int altEnd) {

        ParamUtils.isPositive( startPos, "Genomic positions must be > 0." );
        ParamUtils.isPositive( refEnd, "Genomic positions must be > 0." );
        ParamUtils.isPositive( altEnd, "Genomic positions must be > 0." );

        final int refLength = refEnd - startPos + 1;
        final int altLength = altEnd - startPos + 1;

        // We know it's a frameshift if we have a replacement that is not of a
        // length evenly divisible by 3 because that's how many bases are read at once:
        return ((Math.abs( refLength - altLength ) % 3) != 0);
    }

    /**
     * Determines whether the given reference and alternate alleles constitute an insertion mutation.
     * @param reference The reference {@link Allele}.   Must not be {@code null}.
     * @param alternate The alternate / variant {@link Allele}.   Must not be {@code null}.
     * @return {@code true} if replacing the reference with the alternate results in an insertion.  {@code false} otherwise.
     */
    public static boolean isInsertion(final Allele reference, final Allele alternate) {

        Utils.nonNull(reference);
        Utils.nonNull(alternate);

        // If we have more bases in the alternate, we have an insertion:
        return reference.length() < alternate.length();
    }

    /**
     * Determines whether the given reference and alternate alleles constitute an insertion mutation.
     * @param reference A {@link String} containing the bases of the reference allele.  Must not be {@code null}.
     * @param alternate A {@link String} containing the bases of the alternate / variant allele.  Must not be {@code null}.
     * @return {@code true} if replacing the reference with the alternate results in an insertion.  {@code false} otherwise.
     */
    public static boolean isInsertion(final String reference, final String alternate) {

        Utils.nonNull(reference);
        Utils.nonNull(alternate);

        final String refComparable = reference.replaceAll("\\*", "");
        final String altComperable = alternate.replaceAll("\\*", "");

        // If we have more bases in the alternate, we have an insertion:
        return refComparable.length() < altComperable.length();
    }

    /**
     * Determines whether the given reference and alternate alleles constitute a deletion mutation.
     * @param reference A {@link String} containing the bases of the reference allele.  Must not be {@code null}.
     * @param alternate A {@link String} containing the bases of the alternate / variant allele.  Must not be {@code null}.
     * @return {@code true} if replacing the reference with the alternate results in a deletion.  {@code false} otherwise.
     */
    public static boolean isDeletion(final String reference, final String alternate) {

        Utils.nonNull(reference);
        Utils.nonNull(alternate);

        final String refComparable = reference.replaceAll("\\*", "");
        final String altComperable = alternate.replaceAll("\\*", "");

        // If we have fewer bases in the alternate, we have a deletion:
        return refComparable.length() > altComperable.length();
    }

    /**
     * Determines whether the given reference and alternate alleles constitute a deletion mutation.
     * @param reference The reference {@link Allele}.  Must not be {@code null}.
     * @param alternate The alternate / variant {@link Allele}.  Must not be {@code null}.
     * @return {@code true} if replacing the reference with the alternate results in a deletion.  {@code false} otherwise.
     */
    public static boolean isDeletion(final Allele reference, final Allele alternate) {

        Utils.nonNull(reference);
        Utils.nonNull(alternate);

        // If we have fewer bases in the alternate, we have a deletion:
        return reference.length() > alternate.length();
    }

    /**
     * Determines whether the given reference and alternate alleles constitute an insertion or deletion mutation.
     * @param reference A {@link String} containing the bases of the reference allele.  Must not be {@code null}.
     * @param alternate A {@link String} containing the bases of the alternate / variant allele.  Must not be {@code null}.
     * @return {@code true} if replacing the reference with the alternate results in an insertion or deletion.  {@code false} otherwise.
     */
    public static boolean isIndel(final String reference, final String alternate) {

        Utils.nonNull(reference);
        Utils.nonNull(alternate);

        final String refComparable = reference.replaceAll("\\*", "");
        final String altComperable = alternate.replaceAll("\\*", "");

        // If we do not have the same number of bases in the reference and alternate alleles,
        // then we have an indel:
        return refComparable.length() != altComperable.length();
    }

    /**
     * Determines whether the given reference and alternate alleles constitute an insertion or deletion mutation.
     * @param reference The reference {@link Allele}.  Must not be {@code null}.
     * @param alternate The alternate / variant {@link Allele}.  Must not be {@code null}.
     * @return {@code true} if replacing the reference with the alternate results in an insertion or deletion.  {@code false} otherwise.
     */
    public static boolean isIndel(final Allele reference, final Allele alternate) {

        Utils.nonNull(reference);
        Utils.nonNull(alternate);

        // If we do not have the same number of bases in the reference and alternate alleles,
        // then we have an indel:
        return reference.length() != alternate.length();
    }

    /**
     * Determines whether the given reference and alternate alleles constitute a polymorphism in one or more nucleotides (XNP).
     * @param reference The reference {@link Allele}.  Must not be {@code null}.
     * @param alternate The alternate / variant {@link Allele}.  Must not be {@code null}.
     * @return {@code true} if replacing the reference with the alternate results in an XNP.  {@code false} otherwise.
     */
    public static boolean isXnp(final Allele reference, final Allele alternate) {

        Utils.nonNull(reference);
        Utils.nonNull(alternate);

        // If we have an equal number of bases in the reference and the alternate, we have an ONP:
        return (reference.length() == alternate.length()) && (!reference.equals(alternate));
    }

    /**
     * Determines whether the given reference and alternate alleles constitute a polymorphism in one or more nucleotides (XNP).
     * @param reference A {@link String} containing the bases of the reference allele.  Must not be {@code null}.
     * @param alternate A {@link String} containing the bases of the alternate / variant allele.  Must not be {@code null}.
     * @return {@code true} if replacing the reference with the alternate results in an XNP.  {@code false} otherwise.
     */
    public static boolean isXnp(final String reference, final String alternate) {

        Utils.nonNull(reference);
        Utils.nonNull(alternate);

        final String refComparable = reference.replaceAll("\\*", "");
        final String altComperable = alternate.replaceAll("\\*", "");

        // If we have an equal number of bases in the reference and the alternate, we have an ONP:
        return ((refComparable.length() == altComperable.length()) && (!refComparable.equals(altComperable)));
    }

    /**
     * @param vc {@link VariantContext to test}
     * @return true if the only alternate allele for this VariantContext is a spanning deletion, otherwise false.
     */
    public static boolean isSpanningDeletionOnly(final VariantContext vc){
        return vc.getAlternateAlleles().size() == 1 && GATKVCFConstants.isSpanningDeletion(vc.getAlternateAllele(0));
    }

    /**
     *
     * Attempt to match allele ref/alt pairs, even if the allele pairs in the given variant contexts are equivalent,
     *  but not exact.
     *
     * For example, if variant 1 and variant2 have the same position, but
     * Variant 1: "A", T", "C"
     * Variant 2: "ACCAGGCCCAGCTCATGCTTCTTTGCAGCCTCT", "TCCAGGCCCAGCTCATGCTTCTTTGCAGCCTCT", "A", "AC"
     *
     * Then the returned array would be:  {0, -1}
     * Since A>T matches in variant1 matches the the first alt allele in variant 2.  And A>C does not match any allele
     *  in variant 2.
     *
     * Do not use this method for doing a full split of a variant context into biallelic components.  This method
     *  ignores (and drops) genotypes and INFO attributes.  It is really meant just for alleles, but works on a
     *  VariantContext to better leverage existing code in
     *  {@link GATKVariantContextUtils#trimAlleles(VariantContext, boolean, boolean)}
     *
     * @param variant1 Never {@code null}
     * @param variant2 Never {@code null}
     * @return Returns indexes into the alternate alleles of variant2.  Note that this method assumes that (when biallelic) the variant
     *  contexts are already trimmed.  See {@link GATKVariantContextUtils#trimAlleles(VariantContext, boolean, boolean)}
     *  Never {@code null}.  The length will be equal to the number of alternate alleles in variant1.  A value of -1
     *  indicates no match in variant2.  If the reference alleles do not match, the output array will be populated
     *  exclusively with -1.
     */
    public static int[] matchAllelesOnly(final VariantContext variant1, final VariantContext variant2) {
        Utils.nonNull(variant1);
        Utils.nonNull(variant2);

        // Grab the trivial case:
        if (variant1.isBiallelic() && variant2.isBiallelic()) {
            if (variant1.getAlternateAllele(0).equals(variant2.getAlternateAllele(0)) &&
                    (variant1.getReference().equals(variant2.getReference()))) {
                return new int[]{0};
            } else {
                return new int[]{-1};
            }
        }

        // Handle the case where one or both of the input VCs are not biallelic.
        final int[] result = new int[variant1.getAlternateAlleles().size()];

        // First split (and trim) all variant contexts into biallelics.  We are only going ot be interested in the alleles.
        final List<VariantContext> splitVariants1 = simpleSplitIntoBiallelics(variant1);
        final List<VariantContext> splitVariants2 = simpleSplitIntoBiallelics(variant2);

        // Second, match on ref and alt.  If match occurs add it to the output list.
        for (int i = 0; i < splitVariants1.size(); i++) {
            result[i] = -1;
            for (int j = 0; j < splitVariants2.size(); j++) {
                final VariantContext splitVariant1 = splitVariants1.get(i);
                final VariantContext splitVariant2 = splitVariants2.get(j);
                if (splitVariant1.getAlternateAllele(0).equals(splitVariant2.getAlternateAllele(0))
                        && splitVariant1.getReference().equals(splitVariant2.getReference())) {
                    result[i] = j;
                }
            }
        }

        return result;
    }

    /**
     * Do not use this method for doing a full split of a variant context into biallelic components.  This method
     *  ignores (and drops) genotypes and INFO attributes.  It is really meant just for alleles, but works on a
     *  VariantContext to better leverage existing code in
     *  {@link GATKVariantContextUtils#trimAlleles(VariantContext, boolean, boolean)}
     *
     * This method is trying to be fast, otherwise.
     *
     * @param vc variant context to split into simple biallelics.  Never {@code null}
     * @return a list of variant contexts.  Each will be biallelic.  Length will be the number of alt alleles in the input vc.
     * Note that the variant contexts are usually stripped of attributes and genotypes.  Never {@code null}.  Empty list
     * if variant context has no alt alleles.
     */
    private static List<VariantContext> simpleSplitIntoBiallelics(final VariantContext vc) {
        Utils.nonNull(vc);
        final List<VariantContext> result = new ArrayList<>();

        if (vc.isBiallelic()) {
            return Collections.singletonList(vc);
        } else {
            // Since variant context builders are slow to keep re-creating.  Just create one and spew variant contexts from it, since
            //  the only difference will be the alternate allele.  Initialize the VCB with a dummy alternate allele,
            //  since it will be overwritten in all cases.
            final VariantContextBuilder vcb = new VariantContextBuilder("SimpleSplit", vc.getContig(), vc.getStart(), vc.getEnd(),
                    Arrays.asList(vc.getReference(), Allele.NO_CALL));
            vc.getAlternateAlleles().forEach(allele -> result.add(GATKVariantContextUtils.trimAlleles(
                    vcb.alleles(Arrays.asList(vc.getReference(), allele)).make(true), true, true)
                    )
            );
        }

        return result;
    }

    /**
     * Returns true if the context represents a MNP. If the context supplied contains the GVCF NON_REF symbolic
     * allele, then the determination is performed after discarding the NON_REF symbolic allele.
     *
     * @param vc a Variant context
     * @return true if the context represents an unmixed MNP (i.e. all alleles are non-symbolic and of the same
     *         length > 1), false otherwise
     */
    public static boolean isUnmixedMnpIgnoringNonRef(final VariantContext vc) {
        final List<Allele> alleles = vc.getAlleles();
        int length = vc.getReference().length(); // ref allele cannot be symbolic
        if (length < 2) { return false; }

        for (final Allele a : alleles) {
            if (a.isSymbolic() && !Allele.NON_REF_ALLELE.equals(a)) { return false; }
            else if (!a.isSymbolic() && a.length() != length) { return false; }
        }

        return true;
    }
}

