package bulkprocess

import (
	"archive/zip"
	"bufio"
	"bytes"
	"encoding/csv"
	"fmt"
	"io"
	"log"
	"os"
	"path"
	"strings"

	"github.com/davecgh/go-spew/spew"
	"github.com/jmoiron/sqlx"
)

// Many of the manifest.csv files are invalid. They use commas as delimiters and
// do not use quotes, and yet they have fields (like date) that use commas. Some
// are literally malformed, ending after 1-2 fields, or containing run-ons of
// multiple fields without proper line-breaks. So this requires some processing.
func processInvalidManifest(manifest io.Reader) ([][]string, error) {
	output := [][]string{}

	header := []string{}
	dateField := -1
	i := 0

	// Read over each line of the invalid manifest
	r := bufio.NewReader(manifest)
	for {
		line, err := r.ReadString('\n')
		if err == io.EOF && len(strings.TrimSpace(line)) == 0 {
			// Can't always break if you see an EOF. Specifically, in some
			// files, the final line has no valid ending (neither an EOF nor a
			// newline). In that case, since the entire line still remains
			// despite receiving an EOF signal, so you want to process that.
			// When you then loop around a final time, there will be nothing in
			// the line, just EOF, so it's safe (and necessary) to break.
			break
		} else if err == io.EOF {
			// We hit the end of the file, but we had non-empty data in the
			// line. See if we need to process it. But, after this, we're done.
			// drainBadBuffer = true
		} else if err != nil {
			return nil, err
		}

		// Split into columns, assuming commas are the delimiter
		if i == 0 {
			record := strings.Split(strings.TrimSpace(strings.Replace(line, "discription", "description", -1)), ",")

			header = append(header, record...)
			for i, col := range record {
				if col == "date" {
					dateField = i
					break
				}
			}

			if dateField == -1 {
				return nil, fmt.Errorf("No field named 'date' was found in the manifest")
			}

			output = append(output, header)
			i++

			continue
		}

		// Process non-header records
		record := strings.Split(strings.TrimSpace(line), ",")

		// For now, if there are fewer entries in the field than in the header,
		// we assume it's wrong (Could consider an invalidLayoutX that tries to
		// handle this)
		if len(record) < len(header) {
			log.Println("Skipping record due to being shorter than the header")
			spew.Fdump(os.Stderr, record)
			continue
		}

		// If the invalidLayout1 (basically, assuming a comma in the date field)
		// fails, then we will emit that to stderr and skip the bad record, but
		// continue otherwise.
		parsed, err := invalidLayout1(record, len(header), dateField)
		if err != nil {
			log.Println("Skipping record due to error in parsing:", err.Error())
			spew.Fdump(os.Stderr, record)

			continue
		}

		output = append(output, parsed)

		i++
	}

	return output, nil
}

// invalidLayout1 works for the situation where there are more fields than there
// should be, because there is a comma in the date field.
func invalidLayout1(record []string, headerLen, dateField int) ([]string, error) {
	output := make([]string, headerLen, headerLen)

	for j, col := range record {
		if j == dateField+1 {
			if j < 1 {
				spew.Fdump(os.Stderr, record)
				return output, fmt.Errorf("dateField detected to be %d, nFields %d, not possible to index into %d", dateField, headerLen, j)
			} else if j-1 >= headerLen {
				spew.Fdump(os.Stderr, record)
				return output, fmt.Errorf("dateField detected to be %d, nFields %d, not possible to index into %d", dateField, headerLen, j)
			}
			output[j-1] += ", " + strings.TrimSpace(col)
		} else if j > dateField+1 {
			if j < 1 {
				spew.Fdump(os.Stderr, record)
				return output, fmt.Errorf("dateField detected to be %d, nFields %d, not possible to index into %d", dateField, headerLen, j)
			} else if j-1 >= headerLen {
				spew.Fdump(os.Stderr, record)
				return output, fmt.Errorf("dateField detected to be %d, nFields %d, not possible to index into %d", dateField, headerLen, j)
			}
			output[j-1] = strings.TrimSpace(col)
		} else {
			if j < 0 {
				spew.Fdump(os.Stderr, record)
				return output, fmt.Errorf("dateField detected to be %d, nFields %d, not possible to index into %d", dateField, headerLen, j)
			} else if j >= headerLen {
				spew.Fdump(os.Stderr, record)
				return output, fmt.Errorf("dateField detected to be %d, nFields %d, not possible to index into %d", dateField, headerLen, j)
			}
			output[j] = strings.TrimSpace(col)
		}
	}

	return output, nil
}

func processManifest(manifest io.Reader) ([][]string, error) {
	r := csv.NewReader(manifest)
	return r.ReadAll()
}

func CardiacMRIZipIterator(zipPath string, processOne func(DicomOutput) error) (err error) {
	metadata, err := zipPathToMetadata(zipPath)
	if err != nil {
		return err
	}

	zipName := path.Base(zipPath)

	// First pass: ignore the manifest, look only at the raw dicom metadata.
	// Accumulate all of the metadata we would want so we can use it to update
	// the manifest.

	// Create a lookup map
	dicomMeta := make(map[string]DicomMeta)

	rc, err := zip.OpenReader(zipPath)
	if err != nil {
		return err
	}
	for _, v := range rc.File {
		// Looking only at the dicoms
		if strings.HasPrefix(v.Name, "manifest") {
			continue
		}

		zippedFile, err := v.Open()
		if err != nil {
			return err
		}
		meta, err := DicomToMetadata(zippedFile)
		if err != nil {
			log.Println("Ignoring error and continuing:", err.Error())
			continue
		}

		// Update the lookup
		dicomMeta[v.Name] = *meta
	}

	// Second pass: only look for the manifest.

	// Ropen the file
	rc.Close()
	rc, err = zip.OpenReader(zipPath)
	if err != nil {
		return err
	}
	defer rc.Close()

	for _, v := range rc.File {
		// Looking only at the manifest
		if !strings.HasPrefix(v.Name, "manifest") {
			continue
		}

		zippedFile, err := v.Open()
		if err != nil {
			return err
		}

		// Consume the full manifest into a buffer so we can re-read it if it is
		// invalid
		manifestBuffer := &bytes.Buffer{}
		if _, err := io.Copy(manifestBuffer, zippedFile); err != nil {
			return err
		}

		manifestReader := bytes.NewReader(manifestBuffer.Bytes())

		// Try to process it as a CSV
		fields, err := processManifest(manifestReader)
		if err != nil {

			// If that fails, try again with our manual algorithm
			manifestReader.Seek(0, 0)

			fields, err = processInvalidManifest(manifestReader)
			if err != nil {
				return err
			}
		}

		for loc := range fields {
			if loc == 0 {
				// Discard the header
				continue
			}

			dcm := DicomOutput{}

			dcm, err = stringSliceToDicomStruct(fields[loc])
			if err != nil {
				return err
			}
			dcm.SampleID = metadata.SampleID
			dcm.ZipFile = zipName
			dcm.FieldID = metadata.FieldID
			dcm.Instance = metadata.Instance
			dcm.Index = metadata.Index

			// Fetch the overlay metadata from our lookup table
			if overlaymeta, exists := dicomMeta[dcm.Dicom.Filename]; exists {
				dcm.DicomMeta = overlaymeta
			}

			if dcm.SampleID == "" {
				continue
			}

			if err := processOne(dcm); err != nil {
				return err
			}
		}

		zippedFile.Close()
	}

	return nil
}

// Iterating over the directory will be delegated to the caller.
// Here we will process one zip file.
func ProcessCardiacMRIZip(zipPath string, db *sqlx.DB) (dicoms []DicomOutput, err error) {
	metadata, err := zipPathToMetadata(zipPath)
	if err != nil {
		return nil, err
	}

	zipName := path.Base(zipPath)

	rc, err := zip.OpenReader(zipPath)
	if err != nil {
		return nil, err
	}

	for _, v := range rc.File {
		if !strings.HasPrefix(v.Name, "manifest") {
			continue
		}

		zippedFile, err := v.Open()
		if err != nil {
			return nil, err
		}

		// Consume the full manifest into a buffer so we can re-read it if it is
		// invalid
		manifestBuffer := &bytes.Buffer{}
		if _, err := io.Copy(manifestBuffer, zippedFile); err != nil {
			return nil, err
		}

		manifestReader := bytes.NewReader(manifestBuffer.Bytes())

		// Try to process it as a CSV
		fields, err := processManifest(manifestReader)
		if err != nil {

			// If that fails, try again with our manual algorithm
			manifestReader.Seek(0, 0)

			fields, err = processInvalidManifest(manifestReader)
			if err != nil {
				return nil, err
			}
		}

		dicoms = make([]DicomOutput, len(fields), len(fields))
		for loc := range fields {
			if loc == 0 {
				// Discard the header
				continue
			}

			dicoms[loc], err = stringSliceToDicomStruct(fields[loc])
			if err != nil {
				return nil, err
			}
			dicoms[loc].SampleID = metadata.SampleID
			dicoms[loc].ZipFile = zipName
			dicoms[loc].FieldID = metadata.FieldID
			dicoms[loc].Instance = metadata.Instance
			dicoms[loc].Index = metadata.Index
		}

		zippedFile.Close()
	}

	deleteEmpty := []int{}
	for i := range dicoms {
		if dicoms[i].SampleID == "" {
			deleteEmpty = append(deleteEmpty, i)
		}
	}

	for loc, i := range deleteEmpty {
		dicoms = append(dicoms[:i-loc], dicoms[i+1-loc:]...)
	}

	return dicoms, nil
}

//   Take the manifest and modify its contents to specify the sample,
//     then update the SQLite
//   Annotate with lookup of the axis type *and field ID* if not yet done
//   Take the dicoms and consider renaming them vs keeping the same name.
