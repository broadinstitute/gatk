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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

@CommandLineProgramProperties(
        summary = "Scramble reference panel",
        oneLineSummary = "Scramble reference panel",
        programGroup = VariantManipulationProgramGroup.class
)
@DocumentedFeature
public class ScrambleReferencePanel extends VariantWalker {

    @Argument(shortName = "nSamples",  doc = "number of samples in output", optional = true)
    public int nSamples = -1;

    @Argument(shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    public GATKPath outPath;
    private VariantContextWriter vcfWriter;

    final private LinkedHashMap<String, List<Pair<Integer, Integer>>> samplesToOldIndexMap = new LinkedHashMap<>();

    final private LinkedHashMap<Allele, BitSet> allelesToHaplotypesBitSetMap = new LinkedHashMap<>();
    final private LinkedHashMap<Pair<String, Integer>, BitSet> sharedHaplotypesBitSetMap = new LinkedHashMap<>();

    final private Random rand = new Random();
    private int nextChange;
    private int blockCounter = 0;
    final private int minSharedHaplotypes = 5;

    @Override
    public void onTraversalStart() {
        vcfWriter = createVCFWriter(outPath);

        final VCFHeader inputHeader = getHeaderForVariants();
        final int nSamplesInPanel = inputHeader.getNGenotypeSamples();

        if (nSamples < 0) {
            nSamples = nSamplesInPanel;
        }
        final List<String> sampleNames = new ArrayList<>();
        for (int i=0; i<nSamples; i++) {
            final String sampleName = "sample_" + i;
            samplesToOldIndexMap.put(sampleName, Arrays.asList(Pair.of(i, 0), Pair.of(i, 1)));
            final BitSet bitSet0 = new BitSet(2*nSamplesInPanel);
            bitSet0.set(0, 2*nSamplesInPanel - 1);
            sharedHaplotypesBitSetMap.put(Pair.of(sampleName, 0), bitSet0);
            final BitSet bitSet1 = new BitSet(2*nSamplesInPanel);
            bitSet1.set(0, 2*nSamplesInPanel - 1);
            sharedHaplotypesBitSetMap.put(Pair.of(sampleName, 1), bitSet1);
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

    private void switchHaplotype(final String sampleToChange, final int hapToChange, final int nSamples) {
        final int iSampleToChangeTo = rand.nextInt(nSamples);
        final int hapToChangeTo = rand.nextInt(2);

        final List<Pair<Integer, Integer>> currentHaps = samplesToOldIndexMap.get(sampleToChange);
        currentHaps.set(hapToChange, Pair.of(iSampleToChangeTo, hapToChangeTo));
        //sharedHaplotypesMap.get(Pair.of(sampleToChange, hapToChange)).clear();
        sharedHaplotypesBitSetMap.get(Pair.of(sampleToChange, hapToChange)).set(0, 2*nSamples - 1);
    }

    private boolean checkForRareHaplotype(final String sample, final int iHap, final Allele allele) {
        final BitSet sharedHaplotypesBitSet = sharedHaplotypesBitSetMap.get(Pair.of(sample, iHap));
        final BitSet sharedHaplotypesThisAlleleBitSet = allelesToHaplotypesBitSetMap.get(allele);

        sharedHaplotypesBitSet.and(sharedHaplotypesThisAlleleBitSet);

        return sharedHaplotypesBitSet.cardinality() < minSharedHaplotypes;
    }

    private int getBitSetIndex(final Pair<Integer, Integer> hap) {
        return hap.getLeft() * 2 + hap.getRight();
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext ) {
        final VariantContextBuilder newVC = new VariantContextBuilder(variant.getSource(), variant.getContig(), variant.getStart(), variant.getEnd(), variant.getAlleles());
        if (blockCounter == nextChange) {
            blockCounter = 0;
            nextChange = getNextChange();

            final int iSampleToChange = rand.nextInt(nSamples);
            final int hapToChange = rand.nextInt(2);
            final String sampleToChange = "sample_" + iSampleToChange;
            switchHaplotype(sampleToChange, hapToChange, variant.getNSamples());
        } else {
            blockCounter++;
        }
        final GenotypesContext genotypes = GenotypesContext.create(nSamples);
        final GenotypesContext oldGenotypes = variant.getGenotypes();

        //fill allele to haplotypeMap
//        allelesToHaplotypesMap.clear();
//        for (final Allele allele : variant.getAlleles()) {
//            allelesToHaplotypesMap.put(allele, new HashSet<>());
//        }
//
//        for (int i=0; i<oldGenotypes.size(); i++) {
//            for (int j=0; j<2; j++) {
//                final Allele allele = oldGenotypes.get(i).getAllele(j);
//                allelesToHaplotypesMap.get(allele).add(Pair.of(i,j));
//            }
//        }

        //fill allele to haplotype bitsets
        allelesToHaplotypesBitSetMap.clear();
        for (final Allele allele : variant.getAlleles()) {
            allelesToHaplotypesBitSetMap.put(allele, new BitSet(variant.getNSamples()));
        }

        for (int i=0; i<oldGenotypes.size(); i++) {
            for (int j=0; j<2; j++) {
                final Allele allele = oldGenotypes.get(i).getAllele(j);
                allelesToHaplotypesBitSetMap.get(allele).set(getBitSetIndex(Pair.of(i,j)));
            }
        }

        //check for rare haplotypes
        for (final String sample : samplesToOldIndexMap.keySet()) {
            final List<Pair<Integer, Integer>> currentHaps = samplesToOldIndexMap.get(sample);
            for (int i=0; i<2; i++) {
                Pair<Integer, Integer> hap = currentHaps.get(i);
                Allele allele = oldGenotypes.get(hap.getLeft()).getAllele(hap.getRight());

                while (checkForRareHaplotype(sample, i, allele)) {
                    //haplotype is rare so need to switch
                    switchHaplotype(sample, i, variant.getNSamples());
                    hap = samplesToOldIndexMap.get(sample).get(i);
                    allele = oldGenotypes.get(hap.getLeft()).getAllele(hap.getRight());
                }

            }
        }

        for (final Map.Entry<String, List<Pair<Integer, Integer>>> entry : samplesToOldIndexMap.entrySet()) {
            final List<Allele> alleles = new ArrayList<>();
            for (final Pair<Integer, Integer> hap : entry.getValue()) {
                final Genotype genotype = oldGenotypes.get(hap.getLeft());
                final Allele allele = genotype.getAllele(hap.getRight());
                alleles.add(allele);
            }

            final Genotype newGenotype = new GenotypeBuilder(entry.getKey(),alleles).phased(true).make();
            genotypes.add(newGenotype);
        }
        newVC.genotypes(genotypes);

        vcfWriter.add(newVC.make());
    }

    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }
}


