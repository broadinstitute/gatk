package bulkprocess

import (
	"fmt"
	"image"
	"image/color"
	"io"
	"io/ioutil"
	"log"

	"github.com/gradienthealth/dicom"
	"github.com/gradienthealth/dicom/dicomtag"
)

// DicomMeta holds a small subset of the available metadata which we consider to
// be useful from dicom images.
type DicomOverlayOpts struct {
	TopLeft struct {
		X int
		Y int
	}
	BottomRight struct {
		X int
		Y int
	}
}

func (d DicomOverlayOpts) Subset() bool {
	if d.TopLeft.X == d.BottomRight.X && d.TopLeft.Y == d.BottomRight.Y {
		return false
	}

	return true
}

// DicomToOverlayImage takes in a dicom file (as a reader), a blank image, and
// options, and updates the image according to those options.
func DicomToOverlayImage(dicomReader io.Reader, opts DicomOverlayOpts) ([]image.Image, error) {
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
			{Group: 0x6002, Element: 0x3000},

			// Number of overlay rows
			{Group: 0x6000, Element: 0x0010},

			// Number of overlay columns
			{Group: 0x6000, Element: 0x0011},

			// Number of overlay frames
			{Group: 0x6000, Element: 0x0015},
			{Group: 0x6000, Element: 0x0014},

			// Bits allocated in the overlay
			{Group: 0x6000, Element: 0x0100},
		},
	})
	if parsedData == nil || err != nil {
		return nil, fmt.Errorf("Error reading dicom: %v", err)
	}

	// Determine the number of rows and columns.
	var nCols int
	var nRows int
	for _, elem := range parsedData.Elements {
		if elem.Tag.Compare(dicomtag.Tag{Group: 0x6000, Element: 0x0010}) == 0 {
			nRows = int(elem.Value[0].(uint16))
		}

		if elem.Tag.Compare(dicomtag.Tag{Group: 0x6000, Element: 0x0010}) == 0 {
			nCols = int(elem.Value[0].(uint16))
		}
	}

	if nCols == 0 {
		return nil, fmt.Errorf("Could not determine number of columns")
	}

	var output []image.Image
	img := image.NewGray16(image.Rect(0, 0, nCols, nRows))

	// Determine whether we are using a framed subset of the overlay or the
	// whole thing, and resize our output image if the former.
	needsOffset := false
	if opts.Subset() {
		needsOffset = true

		windowCols := opts.BottomRight.X - opts.TopLeft.X
		windowRows := opts.BottomRight.Y - opts.TopLeft.Y

		img = image.NewGray16(image.Rect(0, 0, windowCols, windowRows))
	}

	// Now iterate over the overlay and generate the image as appropriate.
	for _, elem := range parsedData.Elements {
		if elem.Tag.Compare(dicomtag.Tag{Group: 0x6000, Element: 0x3000}) != 0 {
			continue
		}

		log.Println("Found the Overlay")

		// We're in the overlay data
		for _, enclosed := range elem.Value {
			// There should be one enclosure, and it should contain a slice of
			// bytes, one byte per pixel.

			cellVals, ok := enclosed.([]byte)
			if !ok {
				continue
			}

			n_bits := 8

			// Fill an array with zeroes, sized the nRows * nCols ( == n_bits *
			// len(cellVals) )
			written := make([]int, n_bits*len(cellVals), n_bits*len(cellVals))

			log.Println("Created a", len(written), "array to hold the output")

			for i := range cellVals {
				byte_as_int := cellVals[i]
				for j := 0; j < n_bits; j++ {
					written[i*n_bits+j] = int((byte_as_int >> uint(j)) & 1)
				}
			}

			if needsOffset {
				log.Println("Bounding to", "(", opts.TopLeft.X, opts.TopLeft.Y, ")", "(", opts.BottomRight.X, opts.BottomRight.Y, ")")
			}

			// Iterate over the bytes. There will be 1 value for each cell.
			// So in a 1024x1024 overlay, you will expect 1,048,576 cells.
			for i, overlayValue := range written {
				row := i / nCols
				col := i % nCols

				if overlayValue != 0 {
					if needsOffset {
						if col < opts.TopLeft.X || col > opts.BottomRight.X || row < opts.TopLeft.Y || row > opts.BottomRight.Y {
							// Ignore out of bounds pixels
							continue
						}

						img.SetGray16(col-opts.TopLeft.X, row-opts.TopLeft.Y, color.White)
					} else {
						img.SetGray16(col, row, color.White)
					}
				}
			}
		}
	}

	// buf := new(bytes.Buffer)
	// if err := png.Encode(buf, img); err != nil {
	// 	return nil, err
	// }
	// jpeg.Encode(buf, img, &jpeg.Options{Quality: 100})
	output = append(output, img)

	// log.Printf("This is image %s. Overlay? %t. Overlay width %d x height %d.", instanceNumber, hasOverlay, width, height)
	// log.Printf("Patient orientation:%.2f %.2f %.2f\n", patientX, patientY, patientZ)

	return output, nil
}
