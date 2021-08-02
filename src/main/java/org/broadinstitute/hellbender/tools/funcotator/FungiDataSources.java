package org.broadinstitute.hellbender.tools.funcotator;//package org.broadinstitute.hellbender.tools.funcotator;
//
//public class FungiDataSources {
//    public static final String ashbya_gossypyii = "ashbya_gossypii//";
//    public static final String aspergillus_clavatus = "asperg";
//}


// Java program to read and download
// webpage in html file
import java.io.*;
import java.net.URL;
import java.net.MalformedURLException;

public class FungiDataSources {

    public static void DownloadWebPage(String webpage)
    {
        try {

            // Create URL object
            URL url = new URL(webpage);
            BufferedReader readr =
                    new BufferedReader(new InputStreamReader(url.openStream()));

            // Enter filename in which you want to download
            BufferedWriter writer =
                    new BufferedWriter(new FileWriter("Download2.html"));

            // read each line from stream till end
            String line;
            while ((line = readr.readLine()) != null) {
                writer.write(line);
            }

            readr.close();
            writer.close();
            System.out.println("Successfully Downloaded.");
        }

        // Exceptions
        catch (MalformedURLException mue) {
            System.out.println("Malformed URL Exception raised");
        }
        catch (IOException ie) {
            System.out.println("IOException raised");
        }
    }
    public static void main(String args[])
            throws IOException
    {
        String url = "http://ftp.ensemblgenomes.org/pub/current/bacteria/gtf/bacteria_0_collection/acinetobacter_baumannii_aye_gca_000069245/Acinetobacter_baumannii_aye_gca_000069245.ASM6924v1.51";
        DownloadWebPage(url);
    }
}