package main

import (
	"archive/zip"
	"bytes"
	"encoding/csv"
	"flag"
	"fmt"
	"image/png"
	"io"
	"io/ioutil"
	"log"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/broadinstitute/ml4h/go/bulkprocess"
)

var (
	ZipColumn   = -1
	DicomColumn = -1
)

func main() {
	// Consumes a batch list that contains zipfile names and dicom names.
	// Generally, these should be identifiable from the manifest, and this tool
	// attempts to do that for you.

	// Extracts all requested dicoms to a folder. Optionally, extract the dicom
	// to .jpg in the process (and update the emitted manifest).

	var path, tmpPath, batchPath, delimiter, overlayCutpoints string
	var makeJpeg, makeImageFromOverlay bool

	flag.StringVar(&path, "path", "./", "Path where the UKBB bulk .zip files are being held.")
	flag.StringVar(&tmpPath, "tmp", os.TempDir(), "Path to a temporary directory where dicom and zip files can be stored.")
	flag.StringVar(&batchPath, "batch", "", "File containing a first column with the zip name, and the second column with an inner dicom file name.")
	flag.StringVar(&delimiter, "delimiter", "\t", "Field delimiter for your batch file, if not a tab")
	flag.BoolVar(&makeJpeg, "make-jpeg", false, "Convert files to jpeg instead of dicom?")
	flag.BoolVar(&makeImageFromOverlay, "make-jpeg-from-overlay", false, "Convert overlay data to jpeg instead of dicom?")
	flag.StringVar(&overlayCutpoints, "overlay-cutpoints", "", "Subset the overlay to a rectangle defined by topleft.{X,Y}, bottomright.{X,Y}. Pass 4 comma-sep values: topleft X,topleft Y,bottomright X,bottomright Y")

	flag.Parse()

	overlayOpts := bulkprocess.DicomOverlayOpts{}
	if overlayCutpoints != "" {
		splitVals := strings.Split(overlayCutpoints, ",")
		if len(splitVals) != 4 {
			log.Fatalln("If setting overlay-cutpoints, need to pass 4 values: topleftx,toplefty,bottomrightx,bottomrighty")
		}
		var err error
		if overlayOpts.TopLeft.X, err = strconv.Atoi(splitVals[0]); err != nil {
			log.Fatalln(err)
		}

		if overlayOpts.TopLeft.Y, err = strconv.Atoi(splitVals[1]); err != nil {
			log.Fatalln(err)
		}

		if overlayOpts.BottomRight.X, err = strconv.Atoi(splitVals[2]); err != nil {
			log.Fatalln(err)
		}

		if overlayOpts.BottomRight.Y, err = strconv.Atoi(splitVals[3]); err != nil {
			log.Fatalln(err)
		}
	}

	// Setup the temp folder
	tempSubfolder := filepath.Join(tmpPath, RandOrthoglyphs(15))
	os.Mkdir(tempSubfolder, os.ModePerm)

	log.Println("Copying requested files to a temporary directory:", tempSubfolder)

	// Setup a data structure that enables us to dip into each zip file at most once
	// zipname => [dicomname => struct{}]
	hierarchy := make(map[string]map[string]struct{})

	// Read the batchfile that tells us which samples we're going to keep
	f, err := os.Open(batchPath)
	if err != nil {
		log.Fatalln(err)
	}
	r := csv.NewReader(f)
	r.Comma = []rune(delimiter)[0]

	log.Println("Delimiter is", strconv.QuoteRune(r.Comma))

	var updatedBatchFile *os.File
	if makeJpeg {
		updatedBatchPath := filepath.Join(tmpPath, filepath.Base(batchPath))
		updatedBatchFile, err = os.OpenFile(updatedBatchPath, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0666)
		if err != nil {
			log.Fatalln("Couldn't create the updated batch file:", err.Error())
		}
	}

	// Find all of the files you'd like to keep in the batch
	requestCount := 0
	for {
		record, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatal(err)
		}

		if len(record) < 2 {
			log.Fatalf("Expected at least 2 columns; found %d\n", len(record))
		}

		if requestCount == 0 {
			for i, v := range record {
				if strings.Contains(v, ".zip") {
					ZipColumn = i
					continue
				}
				if strings.Contains(v, ".dcm") {
					DicomColumn = i
					continue
				}
			}

			log.Println("Zero-based Zip field is", ZipColumn, "E.g.:", record[ZipColumn])
			log.Println("Zero-based Dicom field is", DicomColumn, "E.g.:", record[DicomColumn])

			if ZipColumn == -1 || DicomColumn == -1 {
				log.Fatalln("Delimiter could not be detected")
			}
		}

		// Associate the inner dicoms with their outer zips
		dicomMap, exists := hierarchy[record[ZipColumn]]
		if !exists {
			dicomMap = make(map[string]struct{})
		}
		dicomMap[record[DicomColumn]] = struct{}{}
		hierarchy[record[ZipColumn]] = dicomMap

		requestCount++

		if makeJpeg {
			record[DicomColumn] += ".jpg"
			if err := AppendOpenFile(updatedBatchFile, record, r.Comma); err != nil {
				log.Fatalln(err)
			}
		}
	}

	for zf, df := range hierarchy {
		// Dicom list is constructed; now iterate over the zip files

		rc, err := zip.OpenReader(filepath.Join(path, zf))
		if err != nil {
			log.Fatalln(err)
		}

		toFind := len(df)
		for _, v := range rc.File {
			// Iterate over all of the dicoms in the zip
			if toFind == 0 {
				break
			}

			if _, exists := df[v.Name]; !exists {
				continue
			}

			dicomReader, err := v.Open()
			if err != nil {
				log.Fatalln(err)
			}

			outputFileName := v.Name

			if makeJpeg {
				outputFileName += ".jpg"
				jpegs, err := DicomToJpeg(dicomReader)
				if err != nil {
					log.Fatalln(err)
				}

				if len(jpegs) != 1 {
					log.Fatalf("Expected 1 jpeg per dicom, found %d in %s => %s\n", len(jpegs), zf, v.Name)
				}

				// Replace dicomreader so that it now contains our jpeg
				dicomReader.Close()
				dicomReaderPartial := bytes.NewReader(jpegs[0])
				dicomReader = ioutil.NopCloser(dicomReaderPartial)

				// Moved this inside the check. Unclear to me why I had this
				// outside the makeJpg check.
				extractedPath := filepath.Join(tempSubfolder, outputFileName)
				extractedFile, err := os.OpenFile(extractedPath, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, v.Mode())
				if err != nil {
					log.Fatalln(err)
				}

				_, err = io.Copy(extractedFile, dicomReader)
				dicomReader.Close()
			}

			if makeImageFromOverlay {
				outputFileName += ".overlay.png"
				pngs, err := bulkprocess.DicomToOverlayImage(dicomReader, overlayOpts)
				if err != nil {
					log.Fatalln(err)
				}

				if len(pngs) != 1 {
					log.Fatalf("Expected 1 overlay per dicom, found %d in %s => %s\n", len(pngs), zf, v.Name)
				}

				// Replace dicomreader so that it now contains our jpeg
				// dicomReader.Close()
				// dicomReaderPartial := bytes.NewReader(pngs[0])
				// dicomReader = ioutil.NopCloser(dicomReaderPartial)

				extractedPath := filepath.Join(tempSubfolder, outputFileName)
				extractedFile, err := os.OpenFile(extractedPath, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, v.Mode())
				if err != nil {
					log.Fatalln(err)
				}

				err = png.Encode(extractedFile, pngs[0])
				if err != nil {
					log.Fatalln(err)
				}

				dicomReader.Close()

				// By emitting the output folder, we facilitate making this into
				// a pipe-able tool.
				fmt.Println(extractedPath)
			}

			toFind--
		}

		rc.Close()
	}
}

func AppendOpenFile(file *os.File, line []string, delimiter rune) error {
	b := strings.Builder{}

	elems := len(line)
	for i, v := range line {
		b.WriteString(v)

		if i < elems-1 {
			// Prevent naked tab at the end
			b.WriteRune(delimiter)
		}
	}
	b.WriteByte('\n')

	_, err := file.WriteString(b.String())

	return err
}
