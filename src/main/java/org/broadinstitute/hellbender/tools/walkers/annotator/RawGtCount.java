package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotationData;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;

@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Counts of genotypes w.r.t. the reference allele: 0/0, 0/*, */*, i.e. all alts lumped together")
public class RawGtCount implements InfoFieldAnnotation, ReducibleAnnotation {
    private static final String SEPARATOR = ",";

    @Override
    public String getPrimaryRawKey() { return GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY; }

    @Override
    public boolean hasSecondaryRawKeys() {
        return false;
    }

    @Override
    public List<String> getSecondaryRawKeys() {
        return null;
    }

    @Override
    public Map<String, Object> annotateRawData(ReferenceContext ref, VariantContext vc, AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        return null;
    }

    @Override
    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME generics here blow up
    public Map<String, Object> combineRawData(List<Allele> allelesList, List<ReducibleAnnotationData<?>> listOfRawData) {
        ReducibleAnnotationData combinedData = new ReducibleAnnotationData(null);

        for (final ReducibleAnnotationData currentValue : listOfRawData) {
            parseRawDataString(currentValue);
            combineAttributeMap(currentValue, combinedData);
        }
        final Map<String, Object> annotations = new HashMap<>();
        String annotationString = makeRawAnnotationString(allelesList, combinedData.getAttributeMap());
        annotations.put(getPrimaryRawKey(), annotationString);
        return annotations;
    }

    private String makeRawAnnotationString(final List<Allele> vcAlleles, final Map<Allele, List<Integer>> perAlleleData) {
        //TODO: We can't calculate the true hom ref count since there is no annotation on hom ref calls that we are combining
        //TODO: For now it's better not to include the incorrect value of 0.
        return "." + SEPARATOR + perAlleleData.get(Allele.NO_CALL).get(1) + SEPARATOR + perAlleleData.get(Allele.NO_CALL).get(2);
    }

    private void parseRawDataString(ReducibleAnnotationData<List<Integer>> myData) {
        myData.putAttribute(Allele.NO_CALL, parseRawDataString(myData.getRawData()));
    }

    private List<Integer> parseRawDataString(String rawDataString) {
        try {
            final String[] parsed = rawDataString.trim().replaceAll(AnnotationUtils.BRACKET_REGEX, "").split(", *");
            if (parsed.length != 3) {
                throw new UserException.BadInput(String.format("Raw value for %s has %d values, expected 3. Annotation value is %s", GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY, parsed.length, rawDataString));
            }
            final int homRefCount = Integer.parseInt(parsed[0]);
            final int hetCount = Integer.parseInt(parsed[1]);
            final int homVarCount = Integer.parseInt(parsed[2]);
            return Arrays.asList(homRefCount, hetCount, homVarCount);
        } catch (final NumberFormatException e) {
            throw new UserException.BadInput("malformed " + GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY + " annotation: " + rawDataString, e);
        }
    }

    private void combineAttributeMap(ReducibleAnnotationData<List<Integer>> toAdd, ReducibleAnnotationData<List<Integer>> combined) {
        if (combined.getAttribute(Allele.NO_CALL) != null) {
            combined.putAttribute(Allele.NO_CALL, Arrays.asList(combined.getAttribute(Allele.NO_CALL).get(0) + toAdd.getAttribute(Allele.NO_CALL).get(0),
                    combined.getAttribute(Allele.NO_CALL).get(1) + toAdd.getAttribute(Allele.NO_CALL).get(1),
                    combined.getAttribute(Allele.NO_CALL).get(2) + toAdd.getAttribute(Allele.NO_CALL).get(2)));
        } else {
            combined.putAttribute(Allele.NO_CALL, toAdd.getAttribute(Allele.NO_CALL));
        }
    }

    @Override
    public Map<String, Object> finalizeRawData(VariantContext vc, VariantContext originalVC) {
        return null;
    }

    @Override
    public List<String> getKeyNames() {
        return getRawKeyNames();
    }

    @Override
    public List<VCFCompoundHeaderLine> getDescriptions() {
        return Arrays.asList(GATKVCFHeaderLines.getInfoLine(getKeyNames().get(0)));
    }

    @Override
    public List<VCFCompoundHeaderLine> getRawDescriptions() {
        final List<VCFCompoundHeaderLine> lines = new ArrayList<>(1);
        for (final String rawKey : getRawKeyNames()) {
            lines.add(GATKVCFHeaderLines.getInfoLine(rawKey));
        }
        return lines;
    }

    @Override
    public Map<String, Object> annotate(ReferenceContext ref, VariantContext vc, AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        return null;
    }
}
