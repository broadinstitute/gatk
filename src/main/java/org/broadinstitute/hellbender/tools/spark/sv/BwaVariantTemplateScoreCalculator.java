package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMHeaderRecordComparator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.io.FileUtils;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaSparkEngine;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.*;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Created by valentin on 5/19/17.
 */
public class BwaVariantTemplateScoreCalculator implements StructuralVariantTemplateHaplotypeScoreCalculator {

    private final JavaSparkContext ctx;

    public BwaVariantTemplateScoreCalculator(final JavaSparkContext ctx) {
        this.ctx = Utils.nonNull(ctx);
    }

    @Override
    public void calculate(final TemplateHaplotypeScoreTable table) {
        final SAMFileHeader samHeader = new SAMFileHeader();
        samHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        final JavaRDD<GATKRead> unmappedReads = ctx.parallelize(table.templates().stream()
                .flatMap(t -> {
                    final List<Template.Fragment> fragments = t.fragments();
                    return fragments.stream().map(f -> f.toUnmappedRead(samHeader, fragments.size() > 1));

                })
                .collect(Collectors.toList())).cache();
        final Broadcast<Map<String, Template>> templatesByName = ctx.broadcast(table.templates()
            .stream().collect(Collectors.toMap(Template::name, t -> t)));
        for (int i = 0; i < table.numberOfHaplotypes(); i++) {
            final int alleleIndex = i;
            final Haplotype haplotype = table.haplotypes().get(i);
            final SAMSequenceDictionary referenceDictionary = new SAMSequenceDictionary();
            final SAMSequenceRecord contig = new SAMSequenceRecord("1", table.haplotypes().get(i).length());
            referenceDictionary.addSequence(contig);
            final File referenceImage = composeReference(table.haplotypes().get(i));
            final BwaSparkEngine bwa = prepareBwaEngine(referenceImage, samHeader, referenceDictionary);
            final JavaRDD<GATKRead> mappedReads = bwa.align(unmappedReads);
            final JavaPairRDD<Template, Iterable<GATKRead>> templatesAndMappedReads = mappedReads.mapToPair(r -> new Tuple2<String, GATKRead>(r.getName(), r)).groupByKey()
                    .mapToPair(p -> new Tuple2<>(templatesByName.getValue().get(p._1()),  p._2()));
            templatesAndMappedReads.foreach(p ->
                table.set(alleleIndex, table.indexOf(p._1()), calculateLikelihood(haplotype, p._1(), p._2()))
            );
            disposeReference(referenceImage);
        }
    }

    private double calculateLikelihood(final Haplotype haplotype, final Template template, final Iterable<GATKRead> gatkReads) {
        return Double.NaN;
    }

    private void disposeReference(final File referenceImage) {
        referenceImage.delete();
    }

    private BwaSparkEngine prepareBwaEngine(final File referenceImage, final SAMFileHeader header, final SAMSequenceDictionary dict) {
        return new BwaSparkEngine(ctx, referenceImage.getAbsolutePath(), header, dict);
    }

    private File composeReference(final Haplotype haplotype) {
        File dir = null;
        File img = null;
        try {
            dir = File.createTempFile("gatk-sv-bwa-tmp-ref-", ".dir");
            img = File.createTempFile("gatk-sv-bwa-ref-img-", ".fa.img");
            final File fasta = new File(dir, "ref.fa");
            final File fai = new File(dir, "ref.fa.fai");

            final PrintWriter fastaWriter = new PrintWriter(new FileWriter(dir));
            fastaWriter.println(">1");
            fastaWriter.flush();
            final long offset = fasta.length();
            int nextIdx = 0;
            while (nextIdx < haplotype.length()) {
                fastaWriter.println(haplotype.getBaseString(nextIdx, Math.min(haplotype.length(), nextIdx += 60)));
            }
            fastaWriter.close();
            final PrintWriter faiWriter = new PrintWriter(new FileWriter(fai));
            faiWriter.println(String.join("\t", "1", "" + haplotype.length(), "" + offset, "60", "61"));
            faiWriter.close();
            while (true) {
                try {
                    final String bwaIndexCommand = String.format("bwa index %s ", fasta.getAbsoluteFile());
                    // need to add build index to the bwa-jni.
                    final int exitCode = Runtime.getRuntime().exec(bwaIndexCommand).waitFor();
                    if (exitCode != 0) {
                        throw new GATKException(
                                String.format("could not create the haplotype index at '%s' with exit-code '%d' and command '%s'", dir.getAbsolutePath(), exitCode, ""));
                    } else {
                        FileUtils.deleteDirectory(dir);
                        break;
                    }
                } catch (final InterruptedException e) {
                    // nothing to do here, just loop again.
                }
            }
            BwaMemIndex.createIndexImage(fasta.getAbsolutePath(), img.getAbsolutePath());
            FileUtils.deleteDirectory(dir);
            return img;
        } catch (final IOException ex) {
            throw new GATKException("problems when creating the haplotype index image at '" + img + "' from '" + dir + "'");
        }
    }


}
