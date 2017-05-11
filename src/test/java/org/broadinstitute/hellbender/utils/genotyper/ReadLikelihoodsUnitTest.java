package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Test code for {@link ReadLikelihoods}
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class ReadLikelihoodsUnitTest {
    private static final double EPSILON = 1e-6;
    private static final int ODD_READ_START = 101;
    private static final int EVEN_READ_START = 1;

    @Test(dataProvider = "dataSets")
    public void testInstantiationAndQuery(final String[] samples, final Allele[] alleles, final Map<String,List<GATKRead>> reads) {
        final ReadLikelihoods<Allele> result = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);

        Assert.assertEquals(result.numberOfSamples(), samples.length);
        Assert.assertEquals(result.numberOfAlleles(), alleles.length);
        Assert.assertEquals(result.samples().size(), samples.length);
        Assert.assertEquals(result.alleles().size(), alleles.length);

        testSampleQueries(samples, reads, result);
        testAlleleQueries(alleles, result);
        testLikelihoodMatrixQueries(samples, result, null);
    }

    @Test(dataProvider = "dataSets")
    public void testLikelihoodFillingAndQuery(final String[] samples, final Allele[] alleles, final Map<String,List<GATKRead>> reads) {
        final ReadLikelihoods<Allele> result = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);
        final double[][][] likelihoods = fillWithRandomLikelihoods(samples, alleles, result);
        testLikelihoodMatrixQueries(samples, result, likelihoods);
    }

    private double[][][] fillWithRandomLikelihoods(final String[] samples, final Allele[] alleles, final ReadLikelihoods<Allele> result) {
        final Random rnd = Utils.getRandomGenerator();
        final double[][][] likelihoods = new double[samples.length][alleles.length][];
        for (int s = 0; s < likelihoods.length; s++) {
            final LikelihoodMatrix<Allele> sampleLikelihoods = result.sampleMatrix(s);
            for (int a = 0; a < likelihoods[s].length; a++) {
                likelihoods[s][a] = new double[result.sampleReadCount(s)];
                for (int r = 0; r < likelihoods[s][a].length; r++)
                    sampleLikelihoods.set(a,r,likelihoods[s][a][r] = -Math.abs(rnd.nextGaussian()));
            }
        }
        return likelihoods;
    }

    @Test(dataProvider = "dataSets")
    public void testBestAlleles(final String[] samples, final Allele[] alleles, final Map<String,List<GATKRead>> reads) {
        final ReadLikelihoods<Allele> original = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);
        fillWithRandomLikelihoods(samples,alleles,original);
        final int numberOfAlleles = alleles.length;
        for (int s = 0; s < samples.length; s++) {
            final int sampleReadCount = original.sampleReadCount(s);
            final LikelihoodMatrix<Allele> sampleMatrix = original.sampleMatrix(s);
            final double[] bestLkArray = new double[sampleReadCount];
            final int[] bestIndexArray = new int[sampleReadCount];
            final double[] confidenceArray = new double[sampleReadCount];
            for (int r = 0; r < sampleReadCount; r++) {
                int bestindexOfAllele = -1;
                double bestAlleleLk = Double.NEGATIVE_INFINITY;
                double secondBestAlleleLk = Double.NEGATIVE_INFINITY;
                for (int a = 0; a < numberOfAlleles; a++) {
                    final double lk = sampleMatrix.get(a,r);
                    if (lk > bestAlleleLk) {
                        secondBestAlleleLk = bestAlleleLk;
                        bestAlleleLk = lk;
                        bestindexOfAllele = a;
                    } else if (lk > secondBestAlleleLk) {
                        secondBestAlleleLk = lk;
                    }
                }
                bestLkArray[r] = bestAlleleLk;
                confidenceArray[r] = bestAlleleLk - secondBestAlleleLk;
                bestIndexArray[r] = bestindexOfAllele;
            }
            final Collection<ReadLikelihoods<Allele>.BestAllele> bestAlleles = original.bestAlleles();
            for (final ReadLikelihoods<Allele>.BestAllele bestAllele : bestAlleles) {
                final int readIndex = original.readIndex(s,bestAllele.read);
                if (readIndex == -1) continue;
                Assert.assertEquals(bestLkArray[readIndex], bestAllele.likelihood);
                Assert.assertEquals(bestAllele.allele, alleles[bestIndexArray[readIndex]]);
                Assert.assertEquals(bestAllele.confidence, confidenceArray[readIndex], EPSILON);
            }
        }
    }

    @Test(dataProvider = "dataSets")
    public void testBestAlleleMap(final String[] samples, final Allele[] alleles, final Map<String,List<GATKRead>> reads) {
        final ReadLikelihoods<Allele> original = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);
        fillWithRandomLikelihoods(samples,alleles,original);
        final Map<Allele,List<GATKRead>> expected = new LinkedHashMap<>(alleles.length);
        for (final Allele allele : alleles)
            expected.put(allele,new ArrayList<>());

        final int numberOfAlleles = alleles.length;
        for (int s = 0; s < samples.length; s++) {
            final int sampleReadCount = original.sampleReadCount(s);
            final LikelihoodMatrix<Allele> sampleMatrix = original.sampleMatrix(s);
            for (int r = 0; r < sampleReadCount; r++) {
                int bestindexOfAllele = -1;
                double bestAlleleLk = Double.NEGATIVE_INFINITY;
                double secondBestAlleleLk = Double.NEGATIVE_INFINITY;
                for (int a = 0; a < numberOfAlleles; a++) {
                    final double lk = sampleMatrix.get(a,r);
                    if (lk > bestAlleleLk) {
                        secondBestAlleleLk = bestAlleleLk;
                        bestAlleleLk = lk;
                        bestindexOfAllele = a;
                    } else if (lk > secondBestAlleleLk) {
                        secondBestAlleleLk = lk;
                    }
                }
                if ((bestAlleleLk - secondBestAlleleLk) > ReadLikelihoods.BestAllele.INFORMATIVE_THRESHOLD)
                    expected.get(alleles[bestindexOfAllele]).add(sampleMatrix.getRead(r));
            }
        }

        final Map<Allele,List<GATKRead>> actual = original.readsByBestAlleleMap();

        Assert.assertEquals(actual.size(), alleles.length);
        for (final Allele allele : alleles) {
            final List<GATKRead> expectedList = expected.get(allele);
            final List<GATKRead> actualList = actual.get(allele);
            final Set<GATKRead> expectedSet = new LinkedHashSet<>(expectedList);
            final Set<GATKRead> actualSet = new LinkedHashSet<>(actualList);
            Assert.assertEquals(actualSet, expectedSet);
        }
    }

    @Test(dataProvider = "dataSets")
    public void testFilterPoorlyModeledReads(final String[] samples, final Allele[] alleles, final Map<String,List<GATKRead>> reads) {
        final ReadLikelihoods<Allele> original = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);

        for (int s = 0; s < samples.length; s++) {
            final int sampleReadCount = original.sampleReadCount(s);
            for (int r = 0; r < sampleReadCount; r++) {
                if ((r & 1) == 0) continue;
                for (int a = 0; a < alleles.length; a++)
                    original.sampleMatrix(s).set(a,r,-10000);
            }
        }

        final ReadLikelihoods<Allele> result = original.copy();
        result.filterPoorlyModeledReads(2.0);

        for (int s = 0; s < samples.length; s++) {
            final int oldSampleReadCount = original.sampleReadCount(s);
            final int newSampleReadCount = result.sampleReadCount(s);
            Assert.assertEquals(newSampleReadCount, (oldSampleReadCount + 1) / 2);
            final LikelihoodMatrix<Allele> newSampleMatrix = result.sampleMatrix(s);
            final LikelihoodMatrix<Allele> oldSampleMatrix = original.sampleMatrix(s);
            for (int r = 0 ; r < newSampleReadCount; r++) {
                Assert.assertEquals(original.readIndex(s, result.sampleReads(s).get(r)), r * 2);
                for (int a = 0; a < alleles.length; a++) {
                    Assert.assertEquals(newSampleMatrix.get(a, r), oldSampleMatrix.get(a, r * 2));
                }
            }
        }
    }

    @Test(dataProvider = "dataSets")
    public void testFilterReadsToOverlap(final String[] samples, final Allele[] alleles, final Map<String,List<GATKRead>> reads) {
        final ReadLikelihoods<Allele> original = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);
        final SimpleInterval evenReadOverlap = new SimpleInterval(SAM_HEADER.getSequenceDictionary().getSequences().get(0).getSequenceName(), EVEN_READ_START, EVEN_READ_START);
        fillWithRandomLikelihoods(samples,alleles,original);
        final ReadLikelihoods<Allele> result = original.copy();
        result.filterToOnlyOverlappingUnclippedReads(evenReadOverlap);
        final double[][][] newLikelihoods = new double[samples.length][alleles.length][];
        for (int s = 0; s < samples.length ; s++)
            for (int a = 0; a < alleles.length; a++) {
                newLikelihoods[s][a] = new double[(original.sampleReadCount(s) + 1) / 2];
                final LikelihoodMatrix<Allele> sampleMatrix = original.sampleMatrix(s);
                for (int r = 0; r < newLikelihoods[s][a].length; r++) {
                    Assert.assertEquals(result.readIndex(s, sampleMatrix.getRead(r << 1)), r);
                    newLikelihoods[s][a][r] = sampleMatrix.get(a, r << 1);
                }
            }
        testLikelihoodMatrixQueries(samples,result,newLikelihoods);
    }

    @Test(dataProvider = "marginalizationDataSets")
    public void testMarginalizationWithOverlap(final String[] samples, final Allele[] alleles, final Map<String,List<GATKRead>> reads, final Map<Allele,List<Allele>> newToOldAlleleMapping) {
        final ReadLikelihoods<Allele> original = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);
        final GenomeLoc evenReadOverlap = locParser.createGenomeLoc(SAM_HEADER.getSequenceDictionary().getSequences().get(0).getSequenceName(),EVEN_READ_START ,EVEN_READ_START );
        fillWithRandomLikelihoods(samples, alleles, original);
        final ReadLikelihoods<Allele> marginalized = original.marginalize(newToOldAlleleMapping,evenReadOverlap);
        Assert.assertNotNull(marginalized);
        Assert.assertEquals(newToOldAlleleMapping.size(), marginalized.numberOfAlleles());
        for (int a = 0; a < marginalized.numberOfAlleles(); a++) {
            final List<Allele> oldAlleles = newToOldAlleleMapping.get(marginalized.getAllele(a));
            Assert.assertNotNull(oldAlleles);
            for (int s = 0; s < samples.length; s++) {
                final LikelihoodMatrix<Allele> oldSmapleLikelihoods = original.sampleMatrix(s);
                final LikelihoodMatrix<Allele> sampleLikelihoods = marginalized.sampleMatrix(s);
                final int sampleReadCount = sampleLikelihoods.numberOfReads();
                final int oldSampleReadCount = oldSmapleLikelihoods.numberOfReads();
                Assert.assertEquals(sampleReadCount, (oldSampleReadCount + 1) / 2);
                for (int r = 0; r < sampleReadCount; r++) {
                    double oldBestLk = Double.NEGATIVE_INFINITY;
                    for (final Allele oldAllele : oldAlleles) {
                        oldBestLk = Math.max(oldSmapleLikelihoods.get(original.indexOfAllele(oldAllele), r << 1), oldBestLk);
                    }
                    Assert.assertEquals(sampleLikelihoods.get(a, r), oldBestLk);
                }
            }
        }
    }

    @Test(dataProvider = "marginalizationDataSets")
    public void testMarginalization(final String[] samples, final Allele[] alleles, final Map<String,List<GATKRead>> reads, final Map<Allele,List<Allele>> newToOldAlleleMapping) {
        final ReadLikelihoods<Allele> original = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);
        fillWithRandomLikelihoods(samples, alleles, original);
        final ReadLikelihoods<Allele> marginalized = original.marginalize(newToOldAlleleMapping);
        Assert.assertNotNull(marginalized);
        Assert.assertEquals(newToOldAlleleMapping.size(), marginalized.numberOfAlleles());
        for (int a = 0; a < marginalized.numberOfAlleles(); a++) {
            final List<Allele> oldAlleles = newToOldAlleleMapping.get(marginalized.getAllele(a));
            Assert.assertNotNull(oldAlleles);
            for (int s = 0; s < samples.length; s++) {
                final LikelihoodMatrix<Allele> oldSmapleLikelihoods = original.sampleMatrix(s);
                final LikelihoodMatrix<Allele> sampleLikelihoods = marginalized.sampleMatrix(s);
                final int sampleReadCount = sampleLikelihoods.numberOfReads();
                final int oldSampleReadCount = oldSmapleLikelihoods.numberOfReads();
                Assert.assertEquals(oldSampleReadCount, sampleReadCount);
                for (int r = 0; r < sampleReadCount; r++) {
                    double oldBestLk = Double.NEGATIVE_INFINITY;
                    for (final Allele oldAllele : oldAlleles) {
                        oldBestLk = Math.max(oldSmapleLikelihoods.get(original.indexOfAllele(oldAllele), r), oldBestLk);
                    }
                    Assert.assertEquals(sampleLikelihoods.get(a, r), oldBestLk);
                }
            }
        }
    }

    @Test(dataProvider = "dataSets")
    public void testNormalizeBestToZero(final String[] samples, final Allele[] alleles, final Map<String,List<GATKRead>> reads) {
        final ReadLikelihoods<Allele> original = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);
        final double[][][] originalLikelihoods = fillWithRandomLikelihoods(samples,alleles,original);
        final ReadLikelihoods<Allele> result= original.copy();
        result.normalizeLikelihoods(true, Double.NEGATIVE_INFINITY);
        testAlleleQueries(alleles,result);
        final int numberOfAlleles = alleles.length;
        final double[][][] newLikelihoods = new double[originalLikelihoods.length][alleles.length][];
        for (int s = 0; s < samples.length; s++) {
            final int sampleReadCount = original.sampleReadCount(s);
            for (int a = 0; a < numberOfAlleles; a++)
                newLikelihoods[s][a] = new double[sampleReadCount];
            for (int r = 0; r < sampleReadCount; r++) {
                double bestLk = originalLikelihoods[s][0][r];
                for (int a = 1; a < numberOfAlleles; a++) {
                    bestLk = Math.max(bestLk, originalLikelihoods[s][a][r]);
                }
                for (int a = 0; a < numberOfAlleles; a++) {
                    newLikelihoods[s][a][r] = originalLikelihoods[s][a][r] - bestLk;
                }
            }
        }
        testLikelihoodMatrixQueries(samples,result,newLikelihoods);
    }

    @Test(dataProvider = "dataSets")
    public void testNormalizeCapWorstLK(final String[] samples, final Allele[] alleles, final Map<String,List<GATKRead>> reads) {
        final ReadLikelihoods<Allele> original = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);
        final double[][][] originalLikelihoods = fillWithRandomLikelihoods(samples,alleles,original);
        final ReadLikelihoods<Allele> result= original.copy();
        result.normalizeLikelihoods(false, - 0.001);
        testAlleleQueries(alleles,result);
        final int numberOfAlleles = alleles.length;
        final double[][][] newLikelihoods = new double[originalLikelihoods.length][alleles.length][];
        for (int s = 0; s < samples.length; s++) {
            final int sampleReadCount = original.sampleReadCount(s);
            for (int a = 0; a < numberOfAlleles; a++)
                newLikelihoods[s][a] = new double[sampleReadCount];
            for (int r = 0; r < sampleReadCount; r++) {
                double bestAltLk = Double.NEGATIVE_INFINITY;
                for (int a = 0; a < numberOfAlleles; a++) {
                    if (alleles[a].isReference())
                        continue;
                    bestAltLk = Math.max(bestAltLk, originalLikelihoods[s][a][r]);
                }
                if (bestAltLk == Double.NEGATIVE_INFINITY)
                    for (int a = 0; a < numberOfAlleles; a++) {
                        newLikelihoods[s][a][r] = originalLikelihoods[s][a][r];
                    }
                else
                    for (int a = 0; a < numberOfAlleles; a++) {
                        newLikelihoods[s][a][r] = Math.max(originalLikelihoods[s][a][r], bestAltLk - 0.001);
                    }
            }
        }
        testLikelihoodMatrixQueries(samples,result,newLikelihoods);
    }

    @Test(dataProvider = "dataSets")
    public void testNormalizeCapWorstLKAndBestToZero(final String[] samples, final Allele[] alleles, final Map<String,List<GATKRead>> reads) {
        final ReadLikelihoods<Allele> original = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);
        final double[][][] originalLikelihoods = fillWithRandomLikelihoods(samples,alleles,original);
        final ReadLikelihoods<Allele> result= original.copy();
        result.normalizeLikelihoods(true, - 0.001);
        testAlleleQueries(alleles,result);
        final int numberOfAlleles = alleles.length;
        final double[][][] newLikelihoods = new double[originalLikelihoods.length][alleles.length][];
        for (int s = 0; s < samples.length; s++) {
            final int sampleReadCount = original.sampleReadCount(s);
            for (int a = 0; a < numberOfAlleles; a++)
                newLikelihoods[s][a] = new double[sampleReadCount];
            for (int r = 0; r < sampleReadCount; r++) {
                double bestAltLk = Double.NEGATIVE_INFINITY;
                double bestLk = Double.NEGATIVE_INFINITY;
                for (int a = 0; a < numberOfAlleles; a++) {
                    bestLk = Math.max(bestLk, originalLikelihoods[s][a][r]);
                    if (alleles[a].isReference())
                        continue;
                    bestAltLk = Math.max(bestAltLk, originalLikelihoods[s][a][r]);
                }
                if (bestAltLk == Double.NEGATIVE_INFINITY)
                    for (int a = 0; a < numberOfAlleles; a++) {
                        newLikelihoods[s][a][r] = originalLikelihoods[s][a][r] - bestLk;
                    }
                else
                    for (int a = 0; a < numberOfAlleles; a++) {
                        newLikelihoods[s][a][r] = Math.max(originalLikelihoods[s][a][r], bestAltLk - 0.001) - bestLk;
                    }
            }
        }
        testLikelihoodMatrixQueries(samples,result,newLikelihoods);
    }


    @Test(dataProvider = "dataSets")
    public void testAddMissingAlleles(final String[] samples, final Allele[] alleles, final Map<String,List<GATKRead>> reads) {
        final ReadLikelihoods<Allele> original = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);
        final double[][][] originalLikelihoods = fillWithRandomLikelihoods(samples,alleles,original);
        final ReadLikelihoods<Allele> result = original.copy();

        // If all the alleles pass are present in the read-likelihoods collection there is no change.
        result.addMissingAlleles(result.alleles(), Double.NEGATIVE_INFINITY);
        testLikelihoodMatrixQueries(samples,result,originalLikelihoods);

        // If the allele list passed is empty there is no effect.
        result.addMissingAlleles(Collections.emptyList(), Double.NEGATIVE_INFINITY);
        testLikelihoodMatrixQueries(samples,result,originalLikelihoods);

        final Allele newOne;
        final Allele newTwo;
        final Allele newThree;

        // We add a single missing.
        result.addMissingAlleles(Arrays.asList(newOne = Allele.create("ACCCCCAAAATTTAAAGGG".getBytes(), false)),-12345.6);
        Assert.assertEquals(result.numberOfAlleles(), original.numberOfAlleles() + 1);

        // We add too more amongst exisisting alleles:
        result.addMissingAlleles(Arrays.asList(newTwo = Allele.create("ATATATTATATTAATATT".getBytes(), false), result.getAllele(1),
                result.getAllele(0), newThree = Allele.create("TGTGTGTATTG".getBytes(), false), Allele.create("ACCCCCAAAATTTAAAGGG".getBytes(), false)),-6.54321);

        Assert.assertEquals(original.numberOfAlleles() + 3, result.numberOfAlleles());

        final List<Allele> expectedAlleles = new ArrayList<>(original.alleles());
        expectedAlleles.add(newOne); expectedAlleles.add(newTwo); expectedAlleles.add(newThree);

        Assert.assertEquals(result.alleles(), expectedAlleles);

        final double[][][] newLikelihoods = new double[originalLikelihoods.length][][];
        for (int s = 0; s < samples.length; s++) {
            newLikelihoods[s] = Arrays.copyOf(originalLikelihoods[s], originalLikelihoods[s].length + 3);
            final int sampleReadCount = original.sampleReadCount(s);
            final int originalnumberOfAlleles = originalLikelihoods[s].length;
            newLikelihoods[s][originalnumberOfAlleles] = new double[sampleReadCount];
            Arrays.fill(newLikelihoods[s][originalnumberOfAlleles], -12345.6);
            newLikelihoods[s][originalnumberOfAlleles+1] = new double[sampleReadCount];
            Arrays.fill(newLikelihoods[s][originalnumberOfAlleles + 1], -6.54321);
            newLikelihoods[s][originalnumberOfAlleles+2] = new double[sampleReadCount];
            Arrays.fill(newLikelihoods[s][originalnumberOfAlleles + 2], -6.54321);
        }
            testLikelihoodMatrixQueries(samples,result,newLikelihoods);
    }


    @Test(dataProvider = "dataSets")
    public void testAddNonRefAllele(final String[] samples, final Allele[] alleles, final Map<String,List<GATKRead>> reads) {
        final ReadLikelihoods<Allele> original = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);
        final double[][][] originalLikelihoods = fillWithRandomLikelihoods(samples,alleles,original);
        final ReadLikelihoods<Allele> result = original.copy();
        result.addNonReferenceAllele(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
        Assert.assertEquals(result.numberOfAlleles(), original.numberOfAlleles() + 1);
        Assert.assertEquals(result.indexOfAllele(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE), result.numberOfAlleles() - 1);
        final double[][][] newLikelihoods = new double[originalLikelihoods.length][][];
        for (int s = 0; s < samples.length; s++) {
            newLikelihoods[s] = Arrays.copyOf(originalLikelihoods[s], originalLikelihoods[s].length + 1);
            final int sampleReadCount = original.sampleReadCount(s);
            final int ordinarynumberOfAlleles = originalLikelihoods[s].length;
            newLikelihoods[s][ordinarynumberOfAlleles] = new double[sampleReadCount];
            for (int r = 0; r < sampleReadCount; r++) {
                double bestLk = newLikelihoods[s][0][r];
                double secondBestLk = Double.NEGATIVE_INFINITY;
                for (int a = 1; a < ordinarynumberOfAlleles; a++) {
                    final double lk = originalLikelihoods[s][a][r];
                    if (lk > bestLk) {
                        secondBestLk = bestLk;
                        bestLk = lk;
                    } else if (lk > secondBestLk) {
                        secondBestLk = lk;
                    }
                }
                final Median median = new Median();
                final List<Double> qualifylingLikelihoods = new ArrayList<>();
                for (int a = 0; a < ordinarynumberOfAlleles; a++) {
                    if (originalLikelihoods[s][a][r] >= bestLk) continue;
                    qualifylingLikelihoods.add(originalLikelihoods[s][a][r]);
                }
                final double medianLikelihood = median.evaluate(qualifylingLikelihoods.stream().mapToDouble(d -> d).toArray());
                // NaN is returned in cases whether there is no elements in qualifyingLikelihoods.
                // In such case we set the NON-REF likelihood to -Inf.
                final double expectedNonRefLk = !Double.isNaN(medianLikelihood) ? medianLikelihood
                        : ordinarynumberOfAlleles <= 1 ? Double.NaN : bestLk;
                newLikelihoods[s][ordinarynumberOfAlleles][r] = expectedNonRefLk;
            }
        }
        testLikelihoodMatrixQueries(samples,result,newLikelihoods);
    }

    private void testLikelihoodMatrixQueries(final String[] samples, final ReadLikelihoods<Allele> result, final double[][][] likelihoods) {
        for (final String sample : samples) {
            final int indexOfSample = result.indexOfSample(sample);
            final int sampleReadCount = result.sampleReadCount(indexOfSample);
            final int numberOfAlleles = result.numberOfAlleles();
            Assert.assertEquals(result.numberOfAlleles(), numberOfAlleles);
            for (int a = 0; a < numberOfAlleles; a++) {
                Assert.assertEquals(result.sampleReadCount(indexOfSample), sampleReadCount);
                for (int r = 0; r < sampleReadCount; r++) {
                    if (Double.isNaN(result.sampleMatrix(indexOfSample).get(a, r))) {
                        Assert.assertTrue(likelihoods != null && Double.isNaN(likelihoods[indexOfSample][a][r]));
                    } else {
                        Assert.assertEquals(result.sampleMatrix(indexOfSample).get(a, r),
                                likelihoods == null ? 0.0 : likelihoods[indexOfSample][a][r], EPSILON);
                    }
                }
            }
        }
    }

    private void testAlleleQueries(Allele[] alleles, ReadLikelihoods<Allele> result) {
        final Set<Integer> alleleIndices = new LinkedHashSet<>();
        for (final Allele allele : alleles) {
            final int indexOfAllele = result.indexOfAllele(allele);
            Assert.assertTrue(indexOfAllele >= 0);
            Assert.assertFalse(alleleIndices.contains(indexOfAllele));
            alleleIndices.add(indexOfAllele);
            Assert.assertSame(allele, alleles[indexOfAllele]);
        }
    }

    private void testSampleQueries(String[] samples, Map<String, List<GATKRead>> reads, ReadLikelihoods<Allele> result) {
        final Set<Integer> sampleIds = new LinkedHashSet<>(samples.length);
        for (final String sample : samples) {
            final int indexOfSample = result.indexOfSample(sample);
            Assert.assertTrue(indexOfSample >= 0);
            Assert.assertFalse(sampleIds.contains(indexOfSample));
            sampleIds.add(indexOfSample);

            final List<GATKRead> sampleReads = result.sampleReads(indexOfSample);
            final Set<GATKRead> sampleReadsSet = new LinkedHashSet<>(sampleReads);
            final List<GATKRead> expectedSampleReadArray = reads.get(sample);
            final Set<GATKRead> expectedSampleReadsSet = new LinkedHashSet<>(expectedSampleReadArray);
            Assert.assertEquals(sampleReadsSet, expectedSampleReadsSet);

            final int sampleReadCount = sampleReads.size();
            for (int r = 0; r < sampleReadCount; r++) {
                Assert.assertSame(sampleReads.get(r), expectedSampleReadArray.get(r));
                final int readIndex = result.readIndex(indexOfSample, sampleReads.get(r));
                Assert.assertEquals(readIndex,r);
            }
        }
    }

    private String[][] SAMPLE_SETS = {
            {"A","B","C"},
            {"A"},
            {"C","A","D","E","Salsa","Gazpacho"},
    };

    private Allele[][] ALLELE_SETS = {
            {Allele.create("A", true), Allele.create("T"), Allele.create("C")},
            {Allele.create("A", true)},
            {Allele.create("ATTTA"), Allele.create("A", true)},
            {Allele.create("A"), Allele.create("AT", true)},
                    {Allele.create("A", false), Allele.create("AT", false)},
    };

    @DataProvider(name="marginalizationDataSets")
    public Object[][] marginalizationDataSets() {
        try {
            final Random rnd = Utils.getRandomGenerator();
            final Object[][] result = new Object[SAMPLE_SETS.length * ALLELE_SETS.length * ALLELE_SETS.length][];
            int nextIndex = 0;
            for (int s = 0; s < SAMPLE_SETS.length; s++) {
                for (int a = 0; a < ALLELE_SETS.length; a++) {
                    for (int b = 0; b < ALLELE_SETS.length; b++) {
                        if (ALLELE_SETS[b].length < ALLELE_SETS[a].length)
                            result[nextIndex++] = new Object[]{SAMPLE_SETS[s], ALLELE_SETS[a],
                                    dataSetReads(SAMPLE_SETS[s], rnd), randomAlleleMap(ALLELE_SETS[a], ALLELE_SETS[b])
                            };
                    }
                }
            }
            return Arrays.copyOf(result, nextIndex);
        }catch (final Throwable e) {
            throw new RuntimeException(e);
        }
    }

    private Map<Allele,List<Allele>> randomAlleleMap(final Allele[] fromAlleles, final Allele[] toAlleles) {
        final Map<Allele,List<Allele>> result = new LinkedHashMap<>(toAlleles.length);
        for (final Allele toAllele : toAlleles )
            result.put(toAllele, new ArrayList<>(fromAlleles.length));
        final ArrayList<Allele> remaining = new ArrayList<>(Arrays.asList(fromAlleles));
        int nextToIndex = 0;
        final Random rnd = Utils.getRandomGenerator();
        for (int i = 0; i < fromAlleles.length; i++) {
            final int fromindexOfAllele = rnd.nextInt(remaining.size());
            result.get(toAlleles[nextToIndex]).add(remaining.remove(fromindexOfAllele));
            nextToIndex = (nextToIndex + 1) % toAlleles.length;
        }
        return result;
    }


    @DataProvider(name="dataSets")
    public Object[][] dataSets() {
        try {
            final Random rnd = Utils.getRandomGenerator();
            final Object[][] result = new Object[SAMPLE_SETS.length * ALLELE_SETS.length][];
            int nextIndex = 0;
            for (int s = 0; s < SAMPLE_SETS.length; s++)
                for (int a = 0; a < ALLELE_SETS.length; a++) {
                    result[nextIndex++] = new Object[]{SAMPLE_SETS[s], ALLELE_SETS[a],
                            dataSetReads(SAMPLE_SETS[s], rnd)
                    };
                }
            return result;
        }catch (final Throwable e) {
            throw new RuntimeException(e);
        }
    }

    private Map<String,List<GATKRead>> dataSetReads(final String[] samples,
                                                         final Random rnd) {
        final Map<String,List<GATKRead>> result = new LinkedHashMap<>(samples.length);
        for (final String sample : samples) {
            final int readCount = rnd.nextInt(100);
            final List<GATKRead> reads = new ArrayList<>(readCount);
            for (int r = 0; r < readCount; r++) {
                final int alignmentStart = (r & 1) == 0 ? EVEN_READ_START : ODD_READ_START;
                reads.add(ArtificialReadUtils.createArtificialRead(SAM_HEADER,
                        "RRR" + sample + "00" + r, 0, alignmentStart, "AAAAA".getBytes(), new byte[]{30, 30, 30, 30, 30}, "5M"));
            }
            result.put(sample,reads);
        }
        return result;
    }

    @Test(dataProvider="readCountsAndnumberOfAllelesDataSkippingNoAlleleAndWithReference")
    public void testInstantiationAndBasicQueries(final int[] readCounts, final int numberOfAlleles, final boolean hasReference) {
        final SampleList sampleList = sampleList(readCounts);

        final AlleleList<Allele> alleleList = alleleList(numberOfAlleles,hasReference);
        final Map<String,List<GATKRead>> sampleToReads = ReadLikelihoodsUnitTester.sampleToReads(sampleList, readCounts);
        final ReadLikelihoods<Allele> subject = new ReadLikelihoods<>(sampleList,alleleList,sampleToReads);

        AlleleListUnitTester.assertAlleleList(subject, alleleList.asListOfAlleles());
        SampleListUnitTester.assertSampleList(subject,sampleList.asListOfSamples());

        if (hasReference) {
            final int referenceIndex = alleleList.indexOfReference();
            Assert.assertTrue(referenceIndex >= 0);
            Assert.assertEquals(alleleList.indexOfReference(), referenceIndex);
        } else {
            Assert.assertEquals(subject.indexOfReference(), -1);
        }

        testLikelihoodMatrixQueries(alleleList, sampleList, sampleToReads, subject, 0.0);
        testAlleleQueries(alleleList, subject);
        testSampleQueries(sampleList, sampleToReads, subject);
    }

    @Test(dataProvider="readCountsAndnumberOfAllelesDataSkippingNoAlleleAndWithReference")
    public void testInstantiationAndBasicQueries1(final int[] readCounts, final int numberOfAlleles, final boolean hasReference) {
        final SampleList sampleList = sampleList(readCounts);

        final AlleleList<Allele> alleleList = alleleList(numberOfAlleles,hasReference);
        final ReadLikelihoods<Allele> subject = new ReadLikelihoods<>(sampleList,alleleList,Collections.emptyMap());

        final Map<String,List<GATKRead>> sampleToReads = ReadLikelihoodsUnitTester.sampleToReads(sampleList, readCounts);
        final double expectedLik = -0.2;
        subject.addReads(sampleToReads, expectedLik);

        AlleleListUnitTester.assertAlleleList(subject, alleleList.asListOfAlleles());
        SampleListUnitTester.assertSampleList(subject,sampleList.asListOfSamples());

        if (hasReference) {
            final int referenceIndex = alleleList.indexOfReference();
            Assert.assertTrue(referenceIndex >= 0);
            Assert.assertEquals(alleleList.indexOfReference(), referenceIndex);
        } else {
            Assert.assertEquals(subject.indexOfReference(), -1);
        }

        testLikelihoodMatrixQueries(alleleList, sampleList, sampleToReads, subject, expectedLik);
        testAlleleQueries(alleleList, subject);
        testSampleQueries(sampleList, sampleToReads, subject);
    }


    @Test(dataProvider="readCountsAndnumberOfAllelesDataSkippingNoLikelihoodsOrNoAlleleAndWithReference")
    public void testLikelihoodWriting(final int[] readCounts, final int numberOfAlleles, final boolean hasReference) {
        final SampleList sampleList = sampleList(readCounts);

        final AlleleList<Allele> alleleList = alleleList(numberOfAlleles,hasReference);
        final Map<String,List<GATKRead>> sampleToReads = ReadLikelihoodsUnitTester.sampleToReads(sampleList, readCounts);
        final ReadLikelihoods<Allele> subject = new ReadLikelihoods<>(sampleList,alleleList,sampleToReads);

        final int numberOfSamples = readCounts.length;
        int totalLikelihoodsSet = 0;
        int expectedLikelihoodsSet = 0;
        for (int s = 0; s < numberOfSamples; s++) {
            expectedLikelihoodsSet += readCounts[s] * numberOfAlleles;
            final LikelihoodMatrix<Allele> matrix = subject.sampleMatrix(s);
            final int readCount = matrix.numberOfReads();
            for (int a = 0; a < numberOfAlleles; a++)
                for (int r = 0; r < readCount; r++)  {
                    final double likelihood = testLikelihood(s, a, r);
                    Assert.assertNotEquals(likelihood, 0); //Paranoia
                    totalLikelihoodsSet++;
                    matrix.set(a,r,likelihood);
                    Assert.assertEquals(matrix.get(a, r), likelihood);
                }

        }
        Assert.assertEquals(totalLikelihoodsSet, expectedLikelihoodsSet);
    }

    private double testLikelihood(final int indexOfSample, final int indexOfAllele, final int readIndex) {
        return - Math.abs(31 * (indexOfSample + 1) + 101 * indexOfAllele + 1009 * readIndex);
    }

    private final Random rnd = Utils.getRandomGenerator();

    private void testLikelihoodMatrixQueries(final AlleleList<Allele> alleles, final SampleList samples,
                                             final Map<String,List<GATKRead>> sampleToReads, ReadLikelihoods<Allele> result, final double expectedLik) {
        for (final String sample : samples.asListOfSamples()) {
            final int indexOfSample = result.indexOfSample(sample);
            final LikelihoodMatrix<Allele> likelihoodMatrix = result.sampleMatrix(indexOfSample);
            final int sampleReadCount = sampleToReads.get(sample).size();
            final List<GATKRead> reads = sampleToReads.get(sample);
            Assert.assertEquals(likelihoodMatrix.numberOfAlleles(), alleles.numberOfAlleles());
            Assert.assertEquals(likelihoodMatrix.numberOfReads(), sampleReadCount);
            for (int a = 0; a < likelihoodMatrix.numberOfAlleles(); a++) {
                Assert.assertEquals(likelihoodMatrix.getAllele(a), alleles.getAllele(a));
                for (int r = 0; r < sampleReadCount; r++) {
                    Assert.assertEquals(likelihoodMatrix.getRead(r), reads.get(r));
                    Assert.assertEquals(likelihoodMatrix.get(a, r), expectedLik);
                }
            }
        }
    }

    private void testAlleleQueries(final AlleleList<Allele> alleles, ReadLikelihoods<Allele> result) {
        final Set<Integer> alleleIndices = new LinkedHashSet<>();
        for (final Allele allele : alleles.asListOfAlleles()) {
            final int indexOfAllele = result.indexOfAllele(allele);
            Assert.assertTrue(indexOfAllele >= 0);
            Assert.assertFalse(alleleIndices.contains(indexOfAllele));
            alleleIndices.add(indexOfAllele);
            Assert.assertSame(allele, alleles.getAllele(indexOfAllele));
        }
    }

    private void testSampleQueries(final SampleList samples, Map<String, List<GATKRead>> reads,
                                   final ReadLikelihoods<Allele> result) {
        final Set<Integer> sampleIds = new LinkedHashSet<>(samples.numberOfSamples());
        for (final String sample : samples.asListOfSamples()) {
            final int indexOfSample = result.indexOfSample(sample);
            Assert.assertTrue(indexOfSample >= 0);
            Assert.assertFalse(sampleIds.contains(indexOfSample));
            sampleIds.add(indexOfSample);

            final List<GATKRead> sampleReads = result.sampleReads(indexOfSample);
            final Set<GATKRead> sampleReadsSet = new LinkedHashSet<>(sampleReads);
            final List<GATKRead> expectedSampleReadArray = reads.get(sample);
            final Set<GATKRead> expectedSampleReadsSet = new LinkedHashSet<>(expectedSampleReadArray);
            Assert.assertEquals(sampleReadsSet, expectedSampleReadsSet);

            final int sampleReadCount = sampleReads.size();
            for (int r = 0; r < sampleReadCount; r++) {
                Assert.assertSame(sampleReads.get(r), expectedSampleReadArray.get(r));
                final int readIndex = result.readIndex(indexOfSample, sampleReads.get(r));
                Assert.assertEquals(readIndex, r);
            }
        }
    }

    private AlleleList<Allele> alleleList(final int numberOfAlleles, final boolean hasReference) {
        final Allele[] alleles = AlleleListUnitTester.generateRandomUniqueAlleles(numberOfAlleles, 100);
        if (hasReference) {
            final int referenceIndex = rnd.nextInt(numberOfAlleles);
            alleles[referenceIndex] = Allele.create(alleles[referenceIndex].getBases(), true);
        }
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(alleles);
        if (alleleList.numberOfAlleles() != alleles.length) {
            throw new SkipException("repeated alleles, should be infrequent");
        }
        return alleleList;
    }

    private SAMFileHeader SAM_HEADER = ArtificialReadUtils.createArtificialSamHeader(10, 0, 1000);
    final GenomeLocParser locParser = new GenomeLocParser(SAM_HEADER.getSequenceDictionary());


    private int[][] READ_COUNTS = {
            {},
            { 100 },
            { 0 },
            { 0, 0, 0 },
            { 1, 0, 1 },
            { 100, 10 , 100},
            { 1000, 10, 100, 20, 23 }
    };

    private int[] ALLELE_COUNTS = { 0, 1, 2, 3, 10, 20 };

    @DataProvider(name="readCountsAndnumberOfAllelesData")
    public Object[][] readCountsAndnumberOfAllelesData() {
        final Object[][] result = new Object[READ_COUNTS.length * ALLELE_COUNTS.length * 2][];
        int index = 0;
        for (final int[] readCounts : READ_COUNTS)
            for (final int numberOfAlleles : ALLELE_COUNTS) {
                result[index++] = new Object[]{ readCounts, numberOfAlleles, false};
                result[index++] = new Object[]{ readCounts, numberOfAlleles, true};
            }
        return result;
    }

    @DataProvider(name="readCountsAndnumberOfAllelesDataSkippingNoAlleleAndWithReference")
    public Object[][] readCountsAndnumberOfAllelesDataSkippingNoAlleleAndWithReference() {
        final Object[][] raw = readCountsAndnumberOfAllelesData();
        final List<Object[]> result = new ArrayList<>(raw.length);
        for (final Object[] paramSet : raw)
            if (!paramSet[2].equals(true) || !paramSet[1].equals(0))
                result.add(paramSet);
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="readCountsAndnumberOfAllelesDataSkippingNoLikelihoodsOrNoAlleleAndWithReference")
    public Object[][] readCountsAndnumberOfAllelesDataSkippingNoLikelihoodsOrNoAlleleAndWithReference() {
        final Object[][] raw = readCountsAndnumberOfAllelesDataSkippingNoAlleleAndWithReference();
        final List<Object[]> result = new ArrayList<>(raw.length);
        for (final Object[] paramSet : raw) {
            final int[] readCounts = (int[]) paramSet[0];
            final long totalReadCount = MathUtils.sum(readCounts);
            if (totalReadCount > 0)
                result.add(paramSet);
        }
        return result.toArray(new Object[result.size()][]);
    }

    private SampleList sampleList(final int[] readCounts) {
        final List<String> samples = new ArrayList<>(readCounts.length);
        for (int i = 0; i < readCounts.length; i++)
            samples.add("SAMPLE_" + i);
        return new IndexedSampleList(samples);
    }

}
