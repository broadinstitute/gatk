package bulkprocess

import (
	"fmt"
	"io"
	"io/ioutil"
	"strconv"

	"github.com/gradienthealth/dicom"
	"github.com/gradienthealth/dicom/dicomtag"
)

// DicomMeta holds a small subset of the available metadata which we consider to
// be useful from dicom images.
type DicomMeta struct {
	HasOverlay      bool
	OverlayFraction float64
	OverlayRows     int
	OverlayCols     int
	InstanceNumber  string
	PatientX        float64
	PatientY        float64
	PatientZ        float64
}

// Takes in a dicom file (in bytes), emit meta-information
func DicomToMetadata(dicomReader io.Reader) (*DicomMeta, error) {
	dcm, err := ioutil.ReadAll(dicomReader)
	if err != nil {
		return nil, err
	}

	p, err := dicom.NewParserFromBytes(dcm, nil)
	if err != nil {
		return nil, err
	}

	parsedData, err := p.Parse(dicom.ParseOptions{
		DropPixelData: true,
		ReturnTags: []dicomtag.Tag{
			// The overlay data
			{Group: 0x6000, Element: 0x3000},

			// Number of overlay rows
			{Group: 0x6000, Element: 0x0010},

			// Number of overlay columns
			{Group: 0x6000, Element: 0x0011},

			// Number of overlay frames
			{Group: 0x6000, Element: 0x0015},
			{Group: 0x6000, Element: 0x0014},

			// Bits allocated in the overlay
			{Group: 0x6000, Element: 0x0100},

			// Other metadata

			// Patient orientation
			{Group: 0x0020, Element: 0x0020},

			// Image Position
			{Group: 0x0020, Element: 0x0032},

			// Image Orientation
			{Group: 0x0020, Element: 0x0037},

			// Slice location
			{Group: 0x0020, Element: 0x0041},

			// Instance number
			{Group: 0x0020, Element: 0x0013},

			// Patient position
			{Group: 0x0018, Element: 0x5100},
		},
	})
	if parsedData == nil || err != nil {
		return nil, fmt.Errorf("Error reading dicom: %v", err)
	}

	output := &DicomMeta{}
	for _, elem := range parsedData.Elements {
		if elem.Tag.Compare(dicomtag.Tag{Group: 0x6000, Element: 0x3000}) == 0 {
			output.HasOverlay = true

			activeCells := 0
			totalCells := 0

			for _, enclosed := range elem.Value {
				cellVals, ok := enclosed.([]byte)
				if !ok {
					continue
				}

				for _, cell := range cellVals {
					totalCells++
					if cell != 0 {
						activeCells++
					}
				}
			}

			if totalCells > 0 {
				output.OverlayFraction = float64(activeCells) / float64(totalCells)
			}
		}

		if elem.Tag.Compare(dicomtag.Tag{Group: 0x6000, Element: 0x0010}) == 0 {
			output.OverlayRows = int(elem.Value[0].(uint16))
		}

		if elem.Tag.Compare(dicomtag.Tag{Group: 0x6000, Element: 0x0011}) == 0 {
			output.OverlayCols = int(elem.Value[0].(uint16))
		}

		if elem.Tag.Compare(dicomtag.Tag{Group: 0x0020, Element: 0x0013}) == 0 {
			output.InstanceNumber = elem.Value[0].(string)
		}

		if elem.Tag.Compare(dicomtag.Tag{Group: 0x0020, Element: 0x0032}) == 0 {
			output.PatientX, err = strconv.ParseFloat(elem.Value[0].(string), 32)
			if err != nil {
				continue
			}
			output.PatientY, err = strconv.ParseFloat(elem.Value[1].(string), 32)
			if err != nil {
				continue
			}
			output.PatientZ, err = strconv.ParseFloat(elem.Value[2].(string), 32)
			if err != nil {
				continue
			}
		}

		// fmt.Printf("%s: VR %s: %v\n", elem.Tag, elem.VR, elem.String())
	}

	// log.Printf("This is image %s. Overlay? %t. Overlay width %d x height %d.", instanceNumber, hasOverlay, width, height)
	// log.Printf("Patient orientation:%.2f %.2f %.2f\n", patientX, patientY, patientZ)

	return output, nil
}
