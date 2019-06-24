package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.STR;

import it.unimi.dsi.fastutil.ints.Int2ObjectLinkedOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.IntLinkedOpenHashSet;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import htsjdk.variant.variantcontext.Allele;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

/**
 * Represents a collection of alleles on a STR variant site.
 */
public class STRAlleleSet extends AbstractList<STRAllele> implements AlleleList<STRAllele> {

    private final byte[] prefix;

    private final byte[] repeatUnit;

    final int referenceRepeatCount;

    private final int minimumRepeatCount;

    private final int maximumRepeatCount;

    private final Int2ObjectMap<STRAllele> alleleByRepeatCount;

    private final STRAllele[] alleleByIndex;

    /**
     * Compose a {@link STRAlleleSet} from a variant context and reference-context.
     *
     * <p>
     *     The input allele list contains all the alleles including the reference; symbolic alleles are ignored.
     * </p>
     * <p>
     *     The input reference base array is used to determine the number of repeats for the reference allele.
     * </p>
     * <p>
     *     This method will return {@code null} if the variant-context is not an STR, meaning that there is the
     *     list of reference and alternative alleles therein cannot be explained a variation in a repeat unit.
     *     For example is mixture between STR and a SNP or the is an allele with a partial repeat.
     * </p>
     *
     * @param alleles the variant context alleles.
     * @param referenceBases reference bases.
     * @return never {@code null}.
     */
    static STRAlleleSet fromAlleleList(final List<? extends Allele> alleles, final byte[] referenceBases) {
        final Allele refAllele = alleles.stream().filter(Allele::isReference)
                .findFirst().orElseThrow(() -> new IllegalArgumentException("the input list does not contain a reference"));
        final byte[] referenceAlleleBases = refAllele.getBases();
        if (referenceAlleleBases.length > referenceBases.length) {
            throw new IllegalArgumentException("the reference context is shorter than the reference allele");
        }
        for (int i = 0; i < referenceAlleleBases.length; i++) {
            if (referenceBases[i] != referenceAlleleBases[i]) {
                throw new IllegalArgumentException(String.format("the reference context bases don't match the reference allele bases: %s != %s", new String(Arrays.copyOf(referenceBases, referenceAlleleBases.length)), new String(referenceAlleleBases)));
            }
        }
        final List<Allele> strAlleles = new ArrayList<>();
        strAlleles.add(refAllele);

        // We collect alternative vs reference length and base sequence differences:
        final SortedSet<Integer> alternativeLengthDifferences = new IntLinkedOpenHashSet(alleles.size());
        final Set<byte[]> insertBaseSet = new HashSet<>(alleles.size());
        for (final Allele allele : alleles) {
            // Qualifying-alternative alleles are all non-reference non-symbolic alleles:
            if (allele.isReference() || allele.isSymbolic() || allele.isNoCall()) {
                continue;
            }
            final byte[] alleleBases = allele.getBases();
            if (alleleBases.length == referenceAlleleBases.length) {
                return null; // we cannot handle snps and mnps.
            }
            for (int i = 0; i < referenceAlleleBases.length && i < alleleBases.length; i++) {
                if (referenceAlleleBases[i] != alleleBases[i]) {
                    return null; // cannot handle mixed snp within indels.
                }
            }
            final byte[] insertBases = alleleBases.length < referenceAlleleBases.length ?
                    Arrays.copyOfRange(referenceAlleleBases, alleleBases.length, referenceAlleleBases.length)
                    : Arrays.copyOfRange(alleleBases, referenceAlleleBases.length, alleleBases.length);
            insertBaseSet.add(insertBases);
            final int referenceLengthDifference = alleleBases.length - referenceAlleleBases.length;
            alternativeLengthDifferences.add(referenceLengthDifference);
        }
        if (insertBaseSet.size() == 0) {
            return null;
        }
        // Out of the bases differences between alternative and reference alleles we obtain the repeat-unit sequence.
        final byte[] repeatedUnit = composeRepeatedUnit(insertBaseSet);
        if (repeatedUnit == null) { // not a pure STR but a mixed variant.
            return null;
        }

        // Figure out the number of repeats in the reference allele.
        final int referenceUnitCount = countRepeats(referenceBases, 1, repeatedUnit);
        final int[] alternativeReadCounts = alternativeLengthDifferences.stream()
                .mapToInt(a -> a)
                .map(a -> referenceUnitCount + a / repeatedUnit.length)
                .toArray();

        final byte[] prefix = new byte[] { referenceBases[0] };
        return new STRAlleleSet(prefix, repeatedUnit, referenceUnitCount, alternativeReadCounts);
    }

    static STRAlleleSet fromReferenceBases(final byte[] referenceBases, final int unitLength) {
        if (unitLength <= 0) {
            throw new IllegalArgumentException("unit length must be greater than 0: " + unitLength);
        } else if (referenceBases == null) {
            throw new IllegalArgumentException("the reference bases cannot be null");
        } else if (referenceBases.length < 1 + unitLength) {
            throw new IllegalArgumentException("the reference bases passes are too small for such a unit length");
        }
        int referenceRepeatCount = 1;
        for (int offset = 1; offset <= referenceBases.length - unitLength; offset += unitLength) {
            if (!Utils.equalRange(referenceBases, 1, referenceBases, offset, unitLength)) {
                break;
            }
            referenceRepeatCount++;
        }

        return new STRAlleleSet(new byte[] {referenceBases[0]},
                Arrays.copyOfRange(referenceBases, 1, unitLength + 1), referenceRepeatCount);

    }

    STRAlleleSet(final byte[] prefix, final byte[] repeatUnit, final int referenceRepeatCount, final int... alternativeRepeatCounts) {
        if (prefix == null)
        if (alternativeRepeatCounts == null) {
            throw new IllegalArgumentException("the alternative count array cannot be null");
        }
        if (repeatUnit == null) {
            throw new IllegalArgumentException("the repeat unit byte array cannot be null");
        }
        if (referenceRepeatCount < 0) {
            throw new IllegalArgumentException("the reference repeat count cannot be negative");
        }
        this.prefix = prefix == null ? new byte[0] : prefix.clone();
        this.repeatUnit = repeatUnit.clone();
        this.referenceRepeatCount = referenceRepeatCount;
        alleleByIndex = new STRAllele[alternativeRepeatCounts.length + 1];
        alleleByRepeatCount = new Int2ObjectLinkedOpenHashMap<>(alternativeRepeatCounts.length + 1); //new new Int2(alternativeRepeatCounts.length + 1);
        int minimumRepeatCount = referenceRepeatCount;
        int maximumRepeatCount = referenceRepeatCount;
        for (int i = 0; i < alternativeRepeatCounts.length; i++) {
            if (alternativeRepeatCounts[i] == referenceRepeatCount) {
                throw new IllegalArgumentException("no alternative repeat count can be the same as the reference");
            } else if (alternativeRepeatCounts[i] < 0) {
                throw new IllegalArgumentException("no repeat count can be less than 0");
            } else if (alleleByRepeatCount.containsKey(alternativeRepeatCounts[i])) {
                throw new IllegalArgumentException("repeated alternetive allele with count: " + alternativeRepeatCounts[i]);
            } else {
                if (alternativeRepeatCounts[i] < minimumRepeatCount) {
                    minimumRepeatCount = alternativeRepeatCounts[i];
                } else if (alternativeRepeatCounts[i] > maximumRepeatCount) {
                    maximumRepeatCount = alternativeRepeatCounts[i];
                }
            }
        }
        this.minimumRepeatCount = minimumRepeatCount;
        this.maximumRepeatCount = maximumRepeatCount;
        alleleByRepeatCount.put(referenceRepeatCount, alleleByIndex[0] = new STRAllele(this, 0, referenceRepeatCount));
        for (int i = 0; i < alternativeRepeatCounts.length; i++) {
            alleleByRepeatCount.put(alternativeRepeatCounts[i], alleleByIndex[i + 1] =
                    new STRAllele(this, i + 1, alternativeRepeatCounts[i]));
        }
    }

    @Override
    public STRAllele get(final int index) {
        if (index < 0 || index >= alleleByIndex.length)
            throw new IllegalArgumentException("invalid allele index");
        return alleleByIndex[index];
    }

    public STRAllele getByRepeatCount(final int repeatCount) {
        return alleleByRepeatCount.get(repeatCount);
    }

    public int getReferenceRepeatCount() {
        return referenceRepeatCount;
    }

    @Override
    public int size() {
        return alleleByIndex.length;
    }

    public int indexOf(final Allele allele) {
        if (allele instanceof STRAllele) {
            return indexOf((STRAllele) allele);
        } else {
            return Stream.of(alleleByIndex).filter(allele::equals).mapToInt(str -> str.index).findFirst().orElse(-1);
        }
    }

    public int indexOf(final STRAllele allele) {
        if (allele == null) {
            return -1;
        } else if (allele.set == this) {
            return allele.index;
        } else if (!allele.set.compatibleWith(this)) {
            return -1;
        } else if (!alleleByRepeatCount.containsKey(allele.repeatCount)) {
            return -1;
        } else {
            return alleleByRepeatCount.get(allele.repeatCount).index;
        }
    }

    private boolean compatibleWith(final STRAlleleSet other) {
        return Arrays.equals(repeatUnit, other.repeatUnit) && Arrays.equals(prefix, other.prefix);
    }

    byte[] composeRepeatBases(final int repeatCount) {
        final byte[] result = Arrays.copyOf(prefix, repeatUnit.length * (repeatCount - minimumRepeatCount) + prefix.length);
        for (int offset = prefix.length; offset < result.length; offset += repeatUnit.length) {
            System.arraycopy(repeatUnit, 0, result, offset, repeatUnit.length);
        }
        return result;
    }

    /**
     * Returns a human readable representation of the STR allele set.
     * <p>
     *     The format is:
     * </p>
     * <pre>
     *      prefix-bases '(' repeat-unit-bases ')' '{' comma-separated-allele-repeat-count-list '}'
     * </pre>
     * <p>
     *     For example:
     * </p>
     * <pre>
     *     A(TC){2,0,1,4}
     * </pre>
     * <p>
     *     reads as a STR allele where the repeat unit is TC which might be repeated 2, 0 or 4 times. This repeat
     *     is preceded by an A. The first repeat-count (in this case 0) corresponds to the reference allele.
     * </p>
     * <p>
     *     So an equivalent VCF's entry would start like this:
     * </p>
     *
     * <pre>
     *     chr1 555    s555555   ATCTC   A,ATC,ATCTCTCTC
     * </pre>
     *
     * @return never {@code null}.
     */
    @Override
    public String toString() {
        return new String(prefix) + "(" + new String(repeatUnit)
                + "){" + Utils.join(",", Stream.of(alleleByIndex).map(a -> a.repeatCount).toArray()) + "}";
    }

    private static byte[] composeRepeatedUnit(final Set<byte[]> inserts) {
        final Iterator<byte[]> it = inserts.iterator();
        // we find out the shortest insert amongst all:
        byte[] shortestInsert = it.next();
        while (it.hasNext()) {
            final byte[] candidate = it.next();
            if (candidate.length < shortestInsert.length) {
                shortestInsert = candidate;
            }
        }
        // then we try to find the smallest repeat unit length that divides the shortest inserts in equal parts:
        int unitLength;
     unitLengthLoop: for (unitLength = 1; unitLength < shortestInsert.length; unitLength++) {
            if ((shortestInsert.length % unitLength) != 0) {
                continue;
            }
            for (int offset = unitLength; offset < shortestInsert.length; offset += unitLength) {
                if (!Utils.equalRange(shortestInsert, 0, shortestInsert, offset, unitLength)) {
                    continue unitLengthLoop;
                }
            }
            break;
        }

        final byte[] candidate = shortestInsert.length == unitLength ? shortestInsert : Arrays.copyOf(shortestInsert, unitLength);

        // Finally we check that that minimal unit is compatible with all other inserts.
        for (final byte[] insert : inserts) {
            if ((insert.length % unitLength) != 0) {
                return null;
            }
            for (int offset = 0; offset < insert.length; offset += unitLength) {
                if (!Utils.equalRange(candidate, 0, insert, offset, unitLength)) {
                    return null;
                }
            }
        }
        return candidate;
    }

    private static int countRepeats(final byte[] referenceBases, final int from, final byte[] repeatedUnit) {
        final int repeatedUnitLength = repeatedUnit.length;
        final int to = referenceBases.length - repeatedUnitLength;
        int offset;
        for (offset = from; offset <= to; offset += repeatedUnitLength) {
            if (!Utils.equalRange(referenceBases, offset, repeatedUnit, 0, repeatedUnitLength)) {
                return (offset - from) / repeatedUnitLength;
            }
        }
        return (offset - from) / repeatedUnitLength;
    }

    int getRepeatUnitLength() {
        return repeatUnit.length;
    }

    int getMaximumRepeatCount() {
        return maximumRepeatCount;
    }

    private final static Pattern SPEC_PATTERN = Pattern.compile("^\\s*(.*)\\s*\\(\\s*(\\S*)\\s*\\)\\{\\s*(.*?)\\s*\\}\\s*$");

    static STRAlleleSet fromString(final String str) {
        if (str == null) {
            throw new IllegalArgumentException("the input string cannot be null");
        }
        final Matcher matcher = SPEC_PATTERN.matcher(str);
        if (!matcher.matches()) {
            throw new IllegalArgumentException("invalid spec format: " + str);
        }
        final String alleleLengthListString = matcher.group(3);
        final String[] alleleLengthList = alleleLengthListString.split("\\s*,\\s*");
        final int[] alleleLengths = Stream.of(alleleLengthList)
                .mapToInt(s -> {
                    try {
                        final int result = Integer.parseInt(s);
                        if (result < 0) {
                            throw new IllegalArgumentException("invalid allele length: " + s);
                        }
                        return result;
                    } catch (final NumberFormatException ex) {
                        throw new IllegalArgumentException("invalid allele length: " + s, ex);
                    }
                }).toArray();
        if (matcher.group(2).isEmpty()) {
            throw new IllegalArgumentException("the repeat unit cannot have length 0: " + str);
        } else if (alleleLengths.length == 0) {
            throw new IllegalArgumentException("there is no alleles in the spec: " + str);
        } else {
            return new STRAlleleSet(matcher.group(1).getBytes(), matcher.group(2).getBytes(), alleleLengths[0], Arrays.copyOfRange(alleleLengths, 1, alleleLengths.length));
        }
    }

    @Override
    public boolean equals(final Object other) {
        if (other instanceof STRAlleleSet) {
            return equals((STRAlleleSet) other);
        } else {
            return super.equals(other);
        }
    }

    public boolean equals(final STRAlleleSet other) {
        if (other == null) {
            return false;
        } else if (other.size() != this.size()) {
            return false;
        } else if (!Arrays.equals(other.prefix, this.prefix)) {
            return false;
        } else if (!Arrays.equals(other.repeatUnit, this.repeatUnit)) {
            return false;
        } else if (alleleByRepeatCount.keySet().stream().anyMatch(i -> !other.alleleByRepeatCount.containsKey(i))) {
            return false;
        } else {
            return true;
        }
    }

    /**
     * Given a different reference allele returns the allele that would correspond with a particular repeat count.
     * @param reference the reference allele.
     * @param repeatCount the repeat count.
     * @return
     * @throws IllegalArgumentException if the input repeat count is not found in this set.
     */
    public Allele mapAllele(final Allele reference, final int repeatCount) {
        if (reference == null) throw new IllegalArgumentException("the reference must not be null");
        if (!alleleByRepeatCount.containsKey(repeatCount)) throw new IllegalArgumentException("there is no such a allele with the input repeat count");
        final byte[] referenceBases = reference.getBases();
        if ((referenceBases.length - prefix.length) % repeatUnit.length != 0) {
            throw new IllegalArgumentException("the input reference is incompatible with this STR allele set");
        }
        if (!Utils.equalRange(referenceBases, 0, prefix, 0, prefix.length)) {
            throw new IllegalArgumentException("the input reference is incompatible with this STR allele set");
        }
        for (int offset = prefix.length; offset < referenceBases.length; offset += repeatUnit.length) {
            if (!Utils.equalRange(referenceBases, offset, repeatUnit, 0, repeatUnit.length)) {
                throw new IllegalArgumentException("the input reference is compatible with this STR allele set");
            }
        }
        final int referenceLengthDiff = (repeatCount - referenceRepeatCount) * repeatUnit.length;
        if (referenceLengthDiff + referenceBases.length < prefix.length) {
            throw new IllegalArgumentException("the input reference is too short");
        } else if (referenceLengthDiff < 0) {
            return Allele.create(Arrays.copyOfRange(referenceBases, 0, referenceLengthDiff + referenceBases.length), repeatCount == referenceRepeatCount);
        } else {
            final byte[] returnBases = Arrays.copyOf(referenceBases, referenceBases.length + referenceLengthDiff);
            for (int offset = referenceBases.length; offset < returnBases.length; offset += repeatUnit.length) {
                System.arraycopy(repeatUnit, 0, returnBases, offset, repeatUnit.length);
            }
            return Allele.create(returnBases, repeatCount == referenceRepeatCount);
        }
    }

    /**
     * Given a different reference allele returns the allele that would correspond with a particular repeat count.
     * @param reference the reference allele.
     * @return
     * @throws IllegalArgumentException if the input repeat count is not found in this set.
     */
    public int alleleIndex(final Allele reference, final Allele allele) {
        if (reference == null) throw new IllegalArgumentException("the reference must not be null");
        final byte[] referenceBases = reference.getBases();
        if ((referenceBases.length - prefix.length) % repeatUnit.length != 0) {
            return -1; //throw new IllegalArgumentException("the input reference is incompatible with this STR allele set");
        }
        if (!Utils.equalRange(referenceBases, 0, prefix, 0, prefix.length)) {
            return -1; //throw new IllegalArgumentException("the input reference is incompatible with this STR allele set");
        }
        for (int offset = prefix.length; offset < referenceBases.length; offset += repeatUnit.length) {
            if (!Utils.equalRange(referenceBases, offset, repeatUnit, 0, repeatUnit.length)) {
                return -1; //throw new IllegalArgumentException("the input reference is incompatible with this STR allele set");
            }
        }
        final byte[] alleleBases = allele.getBases();
        final int prefixSize = Math.min(alleleBases.length, referenceBases.length);
        for (int i = 0; i < prefixSize; i++) {
            if (alleleBases[i] != alleleBases[i]) {
                return -1; // the input allele together with the referen represent a mix site not a STR.
            }
        }
        final int alleleRefLengthDiff = alleleBases.length - referenceBases.length;
        if (alleleRefLengthDiff % repeatUnit.length != 0) {
            return -1; // the input allele together with the reference represent a mix site not a STR.
        } else {
            if (alleleRefLengthDiff > 0) {
                for (int offset = referenceBases.length; offset < alleleBases.length; offset += repeatUnit.length) {
                    if (!Utils.equalRange(alleleBases, offset, repeatUnit, 0, repeatUnit.length)) {
                        return -1;
                    }
                }
            }
            return alleleByRepeatCount.containsKey(referenceRepeatCount + (alleleRefLengthDiff / repeatUnit.length)) ?
                    alleleByRepeatCount.get(referenceRepeatCount + alleleRefLengthDiff / repeatUnit.length).index : -1;
        }
    }



    public int getMinimumRepeatCount() {
        return minimumRepeatCount;
    }

    @Override
    public int numberOfAlleles() {
        return size();
    }

    @Override
    public int indexOfAllele(final STRAllele allele) {
        return indexOf(allele);
    }

    @Override
    public STRAllele getAllele(int index) {
        return get(index);
    }

    public String getRepeatUnitString() {
        return new String(repeatUnit);
    }

    public STRAllele getReference() {
        return alleleByRepeatCount.get(this.referenceRepeatCount);
    }

    public boolean compatibleAllele(final Allele allele) {
        if (!Utils.equalRange(prefix, 0, allele.getBases(), 0, prefix.length)) {
            return false;
        } else if ((allele.getBases().length - prefix.length) % repeatUnit.length != 0) {
            return false;
        } else {
            for (int offset = prefix.length; offset < allele.getBases().length; offset += repeatUnit.length) {
                if (!Utils.equalRange(repeatUnit, 0, allele.getBases(), offset, repeatUnit.length)) {
                    return false;
                }

            }
            return true;
        }
    }
}

