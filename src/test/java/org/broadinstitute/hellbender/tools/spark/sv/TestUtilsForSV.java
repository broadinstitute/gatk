package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils.getCanonicalChromosomes;

public class TestUtilsForSV {

    public static final ReferenceMultiSource b37_reference = new ReferenceMultiSource(
            GATKBaseTest.b37_reference_20_21, ReferenceWindowFunctions.IDENTITY_FUNCTION);
    public static final SAMSequenceDictionary b37_seqDict = b37_reference.getReferenceSequenceDictionary(null);
    public static final Set<String> b37_canonicalChromosomes = getCanonicalChromosomes(null, b37_seqDict);
    public static final ReferenceMultiSource b38_reference_chr20_chr21 = new ReferenceMultiSource(
            GATKBaseTest.b38_reference_20_21, ReferenceWindowFunctions.IDENTITY_FUNCTION);
    public static final SAMSequenceDictionary b38_seqDict_chr20_chr21 = b38_reference_chr20_chr21.getReferenceSequenceDictionary(null);
    public static final Set<String> b38_canonicalChromosomes = getCanonicalChromosomes(null, b38_seqDict_chr20_chr21);


    // We are having this because SV, especially complex ones, are rare and events on chr20 and 21 are not enough.
    public final static Set<String> hg38FullCanonicalChromosomesSet;
    // this is, as its name suggests, a bare bone seq dict for all chromosomes of hg38: IT CONTAINS ONLY THE CHROMOSOME NAME AND LENGTH
    public final static SAMSequenceDictionary bareBoneHg38SAMSeqDictAllChromosomes;

    static {
        final List<SAMSequenceRecord> hg38Chromosomes = new ArrayList<>();
        final String hg38ChrBareBoneListFile =  GATKBaseTest.toolsTestDir + "/spark/sv/utils/hg38ChrBareBone.txt";
        try (final Stream<String> records = Files.lines(IOUtils.getPath(( hg38ChrBareBoneListFile )))) {
            records.forEach(line -> {
                final String[] fields = line.split("\t", 2);
                hg38Chromosomes.add(new SAMSequenceRecord(fields[0], Integer.valueOf(fields[1])));
            });
            bareBoneHg38SAMSeqDictAllChromosomes = new SAMSequenceDictionary(hg38Chromosomes);
        } catch ( final IOException ioe ) {
            throw new UserException("Can't read nonCanonicalContigNamesFile file " + hg38ChrBareBoneListFile, ioe);
        }

        hg38FullCanonicalChromosomesSet =
                hg38Chromosomes.subList(0, 24).stream().map(SAMSequenceRecord::getSequenceName).collect(Collectors.toCollection(LinkedHashSet::new));
    }
}
