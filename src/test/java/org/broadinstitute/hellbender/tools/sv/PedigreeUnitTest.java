package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.io.File;

public class PedigreeUnitTest extends GATKBaseTest {
    @DataProvider(name = "testTrioProvider")
    public Object[][] testTrioProvider() {
        var relationships = new ArrayList<>();
        for (int i = 0; i < 3; i++) {
            relationships.add(new Pedigree.Relationship(
                    "family_" + i,
                    "child_" + i,
                    "father_" + i,
                    "mother_" + i,
                    (i % 2 == 0) ? Pedigree.Sex.MALE : Pedigree.Sex.FEMALE,
                    ((i + 1) % 3 == 0 ? Pedigree.Phenotype.AFFECTED : Pedigree.Phenotype.UNAFFECTED)
            ));
        }
        return new Object[][] { { relationships } };
    }

    @DataProvider(name = "pedFileProvider")
    public Object[][] pedFileProvider() {
        return new Object[][] {{ publicTestDir + "1kgp.ped" }};
    }

    @Test(dataProvider = "testTrioProvider")
    public void testAddTrioSuccessfully(List<Pedigree.Relationship> relationships) {
        // Arrange
        var pedigree = new Pedigree();

        // Act
        for (var rel : relationships) {
            pedigree.addIndividual(rel);
        }

        // Assert
        Assert.assertEquals(pedigree.getTrioCount(), relationships.size());
        var rel = relationships.get(0);
        Assert.assertEquals(
            pedigree.getTrio(rel.individualId()).toString(),
            String.join("\t", rel.familyId(), rel.individualId(), rel.fatherId(), rel.motherId(), "1", "0"));
    }

    @Test(dataProvider = "pedFileProvider")
    public void testReadAndWriteToPedFileSuccessfully(String pedFilename) throws Exception {
        // Arrange & Act
        var pedigree = new Pedigree(pedFilename);

        // Assert
        Assert.assertEquals(pedigree.getTrioCount(), 6);

        // Arrange
        var tempFile = File.createTempFile("ped_", "_" + System.currentTimeMillis() / 1000L);
        tempFile.deleteOnExit();

        // Act
        pedigree.toPed(tempFile.getAbsolutePath());

        // Assert
        var expLines = Files.readAllLines(Paths.get(pedFilename));
        var outLines = Files.readAllLines(tempFile.toPath());
        Collections.sort(expLines);
        Collections.sort(outLines);
        Assert.assertEquals(outLines.toArray(), expLines.toArray());
    }
}
