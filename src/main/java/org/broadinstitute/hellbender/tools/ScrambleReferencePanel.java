package org.broadinstitute.hellbender.tools;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import picard.cmdline.programgroups.VariantManipulationProgramGroup;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.Set;

@CommandLineProgramProperties(
        summary = "Scramble reference panel",
        oneLineSummary = "Scramble reference panel",
        programGroup = VariantManipulationProgramGroup.class
)
@DocumentedFeature
public class ScrambleReferencePanel extends VariantWalker {

    public int nSamples;

    @Argument(shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    public GATKPath outPath;

    @Argument(shortName = "minShared", optional = true)
    final private int minSharedHaplotypes = 5;

    private VariantContextWriter vcfWriter;

    final private List<Pair<Integer, Integer>> haplotypesBeingCopied = new ArrayList<>();

    private int[] suffixArray;
    private int[] indexNewSubgroups;
    private int[][] sharedHaplotypeIndexBoundaries;
    private int[][] alleleCountsAboveInSuffixArray;
    private int[][] indexNextAlleleOccuranceInSuffixArray;
    final private List<BitSet> alleleHaplotypesBitSets_suffixArrays = new ArrayList<>();
    int currentNAlleles = 0;
    private int[] haplotypeLengths;

    final private Random rand = new Random();
    private int nextChange;
    private int blockCounter = 0;

    @Override
    public void onTraversalStart() {
        vcfWriter = createVCFWriter(outPath);

        final VCFHeader inputHeader = getHeaderForVariants();
        final int nSamplesInPanel = inputHeader.getNGenotypeSamples();
        suffixArray = new int[2*nSamplesInPanel];
        indexNewSubgroups = new int[2*nSamplesInPanel];


        nSamples = nSamplesInPanel;
        final List<String> sampleNames = new ArrayList<>();
        sharedHaplotypeIndexBoundaries = new int[2*nSamplesInPanel][2];
        haplotypeLengths = new int[2*nSamples];
        for (int i=0; i<nSamples; i++) {
            final String sampleName = "sample_" + i;
            for (int j=0; j<2; j++) {
                haplotypesBeingCopied.add(Pair.of(i, j));
                Collections.shuffle(haplotypesBeingCopied);
                sharedHaplotypeIndexBoundaries[2*i+j][0] = 0;
                sharedHaplotypeIndexBoundaries[2*i+j][1] = 2*nSamplesInPanel - 1;
                suffixArray[2*i + j] = 2*i + j;
            }
            sampleNames.add(sampleName);
        }
        final Set<VCFHeaderLine> metaData = inputHeader.getMetaDataInInputOrder();
        final VCFHeader header = new VCFHeader(metaData, sampleNames);
        for (final VCFHeaderLine line : getDefaultToolVCFHeaderLines()) {
            header.addMetaDataLine(line);
        }

        vcfWriter.writeHeader(header);



        if (nSamples > nSamplesInPanel) {
            logger.error("nSamples is greater than the number of samples in the original panel (" + nSamplesInPanel + "), please set to a smaller value");
        }



        nextChange = getNextChange();
    }

    private int getNextChange() {
        return (int) (Math.log(1-rand.nextDouble())/(-0.1));
    }

    private void switchHaplotypeArray(final int hapToChange, final int nSamples) {
        final int hapToChangeTo = rand.nextInt(2*nSamples);
        final Pair<Integer, Integer> hapToSwap1 = haplotypesBeingCopied.get(hapToChange);
        final Pair<Integer, Integer> hapToSwap2 = haplotypesBeingCopied.get(hapToChangeTo);
        haplotypesBeingCopied.set(hapToChange, hapToSwap2);
        haplotypesBeingCopied.set(hapToChangeTo, hapToSwap1);
        sharedHaplotypeIndexBoundaries[hapToChange][0] = 0;
        sharedHaplotypeIndexBoundaries[hapToChange][1] = 2*nSamples - 1;
        haplotypeLengths[hapToChange] = 0;

        sharedHaplotypeIndexBoundaries[hapToChangeTo][0] = 0;
        sharedHaplotypeIndexBoundaries[hapToChangeTo][1] = 2*nSamples - 1;
        haplotypeLengths[hapToChangeTo] = 0;
    }

    private boolean checkForRareHaplotype_suffixArray(final int iHap, final int allele, final int indexNewStart) {
        final int firstIndex = sharedHaplotypeIndexBoundaries[iHap][0];
        final int lastIndex = sharedHaplotypeIndexBoundaries[iHap][1];

//        if(firstIndex<0 || lastIndex<0) {
//            System.out.println("hmm");
//        }
        final int cardinality = alleleCountsAboveInSuffixArray[lastIndex][allele] - alleleCountsAboveInSuffixArray[firstIndex][allele] +
                (indexNextAlleleOccuranceInSuffixArray[lastIndex][allele] == lastIndex? 1 : 0);
        final int newFirstIndexSubgroup = indexNewSubgroups[indexNextAlleleOccuranceInSuffixArray[firstIndex][allele]];
        final int newFirstIndex = newFirstIndexSubgroup >= 0 ? newFirstIndexSubgroup + indexNewStart : -1;
        final int newLastIndex = newFirstIndex >= 0 ? newFirstIndex + cardinality - 1 : -1;
        sharedHaplotypeIndexBoundaries[iHap][0] = newFirstIndex;
        sharedHaplotypeIndexBoundaries[iHap][1] = newLastIndex;

//        if(cardinality < minSharedHaplotypes) {
//            System.out.println("hmm");
//        }

        return cardinality < minSharedHaplotypes;
    }

    private int getBitSetIndex(final int sample, final int chrom) {
        return sample * 2 + chrom;
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext ) {
        final VariantContextBuilder newVC = new VariantContextBuilder(variant.getSource(), variant.getContig(), variant.getStart(), variant.getEnd(), variant.getAlleles());
        if (blockCounter == nextChange) {
            blockCounter = 0;
            nextChange = getNextChange();

            final int iSampleToChange = rand.nextInt(nSamples);
            final int hapToChange = rand.nextInt(2);
            switchHaplotypeArray(getBitSetIndex(iSampleToChange, hapToChange), variant.getNSamples());
        } else {
            blockCounter++;
        }
        final GenotypesContext genotypes = GenotypesContext.create(nSamples);
        final GenotypesContext oldGenotypes = variant.getGenotypes();

//        if (variant.getStart() == 904989) {
//            System.out.println("hmm");
//        }
        //build new suffix array
        alleleHaplotypesBitSets_suffixArrays.clear();
        final List<List<Integer>> newSuffixArray = new ArrayList<>();
        for (int i=0; i<variant.getNAlleles(); i++) {
            newSuffixArray.add(new LinkedList<>());
            alleleHaplotypesBitSets_suffixArrays.add(new BitSet(variant.getNSamples()));
        }



        if (variant.getNAlleles() > currentNAlleles || currentNAlleles > 5) {
            alleleCountsAboveInSuffixArray = new int[2 * variant.getNSamples()][variant.getNAlleles()];
            indexNextAlleleOccuranceInSuffixArray = new int[2 * variant.getNSamples()][variant.getNAlleles()];
            currentNAlleles = variant.getNAlleles();
        }

        final int[] unfilledNextAllelesIndecies = new int[variant.getNAlleles()];
        int prevAlleleIndex = 0;
        for (int i=0; i<suffixArray.length; i++) {
            final int iHap = suffixArray[i];
            final int iSample = iHap/2;
            final int iChrom = iHap%2;
            final Allele allele = oldGenotypes.get(iSample).getAllele(iChrom);
            final int alleleIndex = variant.getAlleleIndex(allele);
            final List<Integer> newSubgroup = newSuffixArray.get(alleleIndex);
            indexNewSubgroups[i] = newSubgroup.size();
            newSubgroup.add(iSample*2 + iChrom);
            alleleHaplotypesBitSets_suffixArrays.get(alleleIndex).set(i);
            if (i>0) {
                alleleCountsAboveInSuffixArray[i]=Arrays.copyOf(alleleCountsAboveInSuffixArray[i-1],variant.getNAlleles());
                alleleCountsAboveInSuffixArray[i][prevAlleleIndex]++;
            }
            for (int j=unfilledNextAllelesIndecies[alleleIndex]; j<=i; j++) {
                indexNextAlleleOccuranceInSuffixArray[j][alleleIndex] = i;
            }
            unfilledNextAllelesIndecies[alleleIndex] = i+1;
            prevAlleleIndex = alleleIndex;
        }

        //fill in end of indexNextAlleleOccuranceInSuffixArray
        for (int i=0; i<variant.getNAlleles(); i++) {
            for (int j=unfilledNextAllelesIndecies[i]; j<indexNextAlleleOccuranceInSuffixArray.length; j++) {
                indexNextAlleleOccuranceInSuffixArray[j][i] = -1;
            }
        }

        //check that there is more than one allele with AC>minSharedHaplotypes
        int allowedAllele = -1;
        int allowedAlleleCount = 0;
        for (int iAllele = 0; iAllele < variant.getNAlleles(); iAllele++) {
            if (alleleCountsAboveInSuffixArray[alleleCountsAboveInSuffixArray.length - 1][iAllele] +
                    ((indexNextAlleleOccuranceInSuffixArray[alleleCountsAboveInSuffixArray.length - 1][iAllele] == alleleCountsAboveInSuffixArray.length - 1) ? 1 : 0) >= minSharedHaplotypes) {
                allowedAllele = iAllele;
                allowedAlleleCount++;
            }
        }

        if (allowedAlleleCount < 2) {
            final Allele onlyAllowedAllele = variant.getAlleles().get(allowedAllele);
            final List<Allele> alleles = Arrays.asList(onlyAllowedAllele, onlyAllowedAllele);
            for(int i=0; i<nSamples; i++) {
                final Genotype newGenotype = new GenotypeBuilder("sample_" + i, alleles).phased(true).make();
                genotypes.add(newGenotype);
            }

            newVC.genotypes(genotypes);
            vcfWriter.add(newVC.make());

            for (int i=0;i<haplotypeLengths.length; i++) {
                haplotypeLengths[i]++;
            }

            return;

        }


        final List<Integer> newSuffixArrayStartIndecies = new ArrayList<>();
        newSuffixArrayStartIndecies.add(0);
        for (int i=1; i<newSuffixArray.size(); i++) {
            newSuffixArrayStartIndecies.add(newSuffixArrayStartIndecies.get(i-1) + newSuffixArray.get(i-1).size());
        }

        for (int iHap=0; iHap<haplotypesBeingCopied.size(); iHap++) {
            Pair<Integer, Integer> currentHap = haplotypesBeingCopied.get(iHap);
            Allele allele = oldGenotypes.get(currentHap.getLeft()).getAllele(currentHap.getRight());
            int alleleIndex = variant.getAlleleIndex(allele);

            while(checkForRareHaplotype_suffixArray(iHap, alleleIndex, newSuffixArrayStartIndecies.get(alleleIndex))) {
                switchHaplotypeArray(iHap, variant.getNSamples());
                currentHap = haplotypesBeingCopied.get(iHap);
                allele = oldGenotypes.get(currentHap.getLeft()).getAllele(currentHap.getRight());
                alleleIndex = variant.getAlleleIndex(allele);
                blockCounter = 0;
            }
        }

        //update suffix array
        int index = 0;
        for(final List<Integer> suffixArrayGroup : newSuffixArray) {
            for (final Integer iHap : suffixArrayGroup) {
                suffixArray[index] = iHap;
                index++;
            }
        }



        for(int i=0; i<nSamples; i++) {
            final List<Allele> alleles = new ArrayList<>();
            for (int j=0; j<2; j++) {
                final int iHap = getBitSetIndex(i, j);
                final Pair<Integer, Integer> hapToCopy = haplotypesBeingCopied.get(iHap);
                final Genotype genotype = oldGenotypes.get(hapToCopy.getLeft());
                final Allele allele = genotype.getAllele(hapToCopy.getRight());
                alleles.add(allele);
            }

            final Genotype newGenotype = new GenotypeBuilder("sample_"+i,alleles).phased(true).make();
            genotypes.add(newGenotype);
        }

        newVC.genotypes(genotypes);

        vcfWriter.add(newVC.make());

        for (int i=0;i<haplotypeLengths.length; i++) {
            haplotypeLengths[i]++;
        }
    }

    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }
}


