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
import java.util.BitSet;
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

    @Argument(shortName = "nSamples",  doc = "number of samples in output", optional = true)
    public int nSamples = -1;

    @Argument(shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    public GATKPath outPath;
    private VariantContextWriter vcfWriter;

    final private List<Pair<Integer, Integer>> haplotypesBeingCopied = new ArrayList<>();

    final private List<BitSet> alleleHaplotypesBitSets = new ArrayList<>();
    final private List<BitSet> sharedHaplotypesBitSets = new ArrayList<>();

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
            for (int j=0; j<2; j++) {
                final int iSampleToCopy = rand.nextInt(nSamples);
                final int iChromToCopy = rand.nextInt(2);
                haplotypesBeingCopied.add(Pair.of(iSampleToCopy, iChromToCopy));
                final BitSet bitSet = new BitSet(2*nSamplesInPanel);
                bitSet.set(0, 2*nSamplesInPanel - 1);
                sharedHaplotypesBitSets.add(bitSet);
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
        final int iSampleToChangeTo = rand.nextInt(nSamples);
        final int hapToChangeTo = rand.nextInt(2);

        haplotypesBeingCopied.set(hapToChange, Pair.of(iSampleToChangeTo, hapToChangeTo));
        sharedHaplotypesBitSets.get(hapToChange).set(0, 2*nSamples - 1);
    }

    private boolean checkForRareHaplotype(final int iHap, final int allele) {
        final BitSet sharedHaplotypesBitSet = sharedHaplotypesBitSets.get(iHap);
        final BitSet sharedHaplotypesThisAlleleBitSet = alleleHaplotypesBitSets.get(allele);

        sharedHaplotypesBitSet.and(sharedHaplotypesThisAlleleBitSet);

        return sharedHaplotypesBitSet.cardinality() < minSharedHaplotypes;
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

        //fill allele bitSets
        alleleHaplotypesBitSets.clear();
        for (int i=0; i<variant.getNAlleles(); i++) {
            alleleHaplotypesBitSets.add(new BitSet(variant.getNSamples()));
        }

        for (int i=0; i<oldGenotypes.size(); i++) {
            for (int j=0; j<2; j++) {
                final Allele allele = oldGenotypes.get(i).getAllele(j);
                alleleHaplotypesBitSets.get(variant.getAlleleIndex(allele)).set(getBitSetIndex(i, j));
            }
        }

        //check for rare haplotypes
        for (int iHap=0; iHap<haplotypesBeingCopied.size(); iHap++) {
            Pair<Integer, Integer> currentHap = haplotypesBeingCopied.get(iHap);
            Allele allele = oldGenotypes.get(currentHap.getLeft()).getAllele(currentHap.getRight());

            while(checkForRareHaplotype(iHap, variant.getAlleleIndex(allele))) {
                switchHaplotypeArray(iHap, variant.getNSamples());
                currentHap = haplotypesBeingCopied.get(iHap);
                allele = oldGenotypes.get(currentHap.getLeft()).getAllele(currentHap.getRight());
                blockCounter = 0;
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
    }

    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }
}


