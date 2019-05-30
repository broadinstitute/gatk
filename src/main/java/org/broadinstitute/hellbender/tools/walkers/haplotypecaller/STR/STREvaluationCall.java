package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.STR;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by valentin on 12/4/16.
 */
class STREvaluationCall {
    public Allele reference;
    public List<Allele> alleles;
    public boolean isCalled;
    public double confidence;

    public boolean isHomRef() {
        return alleles.size() > 0 && alleles.stream().allMatch(Allele::isReference);
    }

    public boolean containsIndel() {
        final int refLength = reference.length();
        return alleles.stream().anyMatch(a -> a.length() != refLength);
    }

    public boolean isHet() {
        if (alleles.size() < 2) {
            return false;
        } else if (!alleles.stream().anyMatch(Allele::isReference)) {
            return false;
        } else if (!alleles.stream().anyMatch(Allele::isNonReference)) {
            return false;
        } else {
            return true;
        }
    }

    public Allele hetAlternative() {
        return alleles.stream().filter(Allele::isNonReference).findAny().get();
    }

    public boolean isHomVar() {
        final List<Allele> distinct = alleles.stream().distinct().collect(Collectors.toList());
        return distinct.size() == 1 && !distinct.get(0).isReference();
    }

    public Allele homAlternative() {
        final List<Allele> distinct = alleles.stream().distinct().collect(Collectors.toList());
        return distinct.size() != 1 || distinct.get(0).isReference() ? null : distinct.get(0);
    }

    public boolean isCompositeVar() {
        final List<Allele> distinct = alleles.stream().distinct().collect(Collectors.toList());
        return distinct.size() > 2 && alleles.stream().noneMatch(Allele::isReference);
    }

    public Collection<Object> alternatives() {
        return alleles.stream().filter(Allele::isNonReference).distinct().collect(Collectors.toList());
    }

    static Pair<STREvaluationCall, STREvaluationCall> harmonize(final STREvaluationCall left, final STREvaluationCall right) {
        if (left.reference.equals(right.reference)) {
            return new ImmutablePair<>(left, right);
        } else if (left.reference.length() > right.reference.length()) {
            final byte[] suffix = Arrays.copyOfRange(left.reference.getBases(), right.reference.getBases().length, left.reference.getBases().length);
            final List<Allele> otherAlleles = right.alleles.stream()
                    .map(a -> Allele.create(Utils.concat(a.getBases(), suffix), a.isReference()))
                    .collect(Collectors.toList());
            return new ImmutablePair<>(left, new STREvaluationCall(left.reference, otherAlleles, right.isCalled, right.confidence));
        } else {
            final byte[] suffix = Arrays.copyOfRange(right.reference.getBases(), left.reference.getBases().length, right.reference.getBases().length);
            final List<Allele> otherAlleles = left.alleles.stream()
                    .map(a -> Allele.create(Utils.concat(a.getBases(), suffix), a.isReference()))
                    .collect(Collectors.toList());
            return new ImmutablePair<>(new STREvaluationCall(right.reference, otherAlleles, left.isCalled, left.confidence), right);
        }
    }

    private STREvaluationCall(final Allele reference, final List<Allele> alleles, final boolean isCalled, final double confidence) {
        this.reference = reference;
        this.alleles = alleles;
        this.isCalled = isCalled;
        this.confidence = confidence;
    }

    static STREvaluationCall fromVariantContext(final VariantContext vc, final ReferenceContext ref, final boolean missimgMeansHomRef) {
        if (vc == null) {
            if (missimgMeansHomRef) {
                final Allele refAllele = Allele.create(ref.getBase(), true);
                return new STREvaluationCall(refAllele, Arrays.asList(refAllele, refAllele), true, 99);
            } else {
                return new STREvaluationCall(Allele.create(ref.getBase(), true), Collections.emptyList(), false, 0);
            }
        } else if (vc.getStart() == ref.getInterval().getStart()) {
            Allele reference = vc.getReference();
            List<Allele> alleles = vc.getGenotype(0).getAlleles().stream()
                    .filter(a -> a.isCalled() && !a.isSymbolic())
                    .collect(Collectors.toList());
            byte[] commonSuffix = reference.getBases();
            for (final Allele allele : alleles) {
                int i = allele.length() - 1;
                int j = commonSuffix.length - 1;
                while (i > 0 && j >= 0 && allele.getBases()[i] == commonSuffix[j]) {
                    i--;
                    j--;
                }
                if (j >= 0) {
                    commonSuffix = Arrays.copyOfRange(commonSuffix, j + 1, commonSuffix.length);
                }
            }
            if (commonSuffix.length == reference.length()) {
                reference = Allele.create(reference.getBases()[0], true);
                alleles = Collections.nCopies(2, reference);
            } else if (commonSuffix.length > 0) {
                final int commonSuffixLength = commonSuffix.length;
                reference = Allele.create(Arrays.copyOfRange(reference.getBases(), 0, reference.length() - commonSuffixLength), true);
                alleles = alleles.stream()
                        .map(a -> Allele.create(Arrays.copyOfRange(a.getBases(), 0, a.length() - commonSuffixLength), a.isReference()))
                        .collect(Collectors.toList());
            }
            return new STREvaluationCall(reference, alleles, !alleles.isEmpty(), vc.getPhredScaledQual() < 0 ? (vc.getGenotype(0).hasGQ() ? vc.getGenotype(0).getGQ() :  0) : vc.getPhredScaledQual());
        } else if (AssessSTRs.isNonVariantBlock(vc)) {
            final Allele refAllele = Allele.create(ref.getBase(), true);
            return new STREvaluationCall(refAllele, Arrays.asList(refAllele, refAllele), true, vc.getPhredScaledQual() < 0 ? (vc.getGenotype(0).hasGQ() ? vc.getGenotype(0).getGQ() :  0) : vc.getPhredScaledQual());
        } else {
            return null;
        }
    }
}
