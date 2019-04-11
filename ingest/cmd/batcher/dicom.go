package main

import (
	"bytes"
	"fmt"
	"image"
	"image/color"
	"image/jpeg"
	"io"
	"io/ioutil"

	"github.com/gradienthealth/dicom"
	"github.com/gradienthealth/dicom/dicomtag"
)

// Takes in a dicom file (in bytes), outputs one or more jpeg file equivalents
// (in bytes)
func DicomToJpeg(dicomReader io.Reader) ([][]byte, error) {
	dcm, err := ioutil.ReadAll(dicomReader)
	if err != nil {
		return nil, err
	}

	p, err := dicom.NewParserFromBytes(dcm, nil)
	if err != nil {
		return nil, err
	}

	parsedData, err := p.Parse(dicom.ParseOptions{DropPixelData: false})
	if parsedData == nil || err != nil {
		return nil, fmt.Errorf("Error reading dicom: %v", err)
	}

	var output [][]byte

	for _, elem := range parsedData.Elements {
		if elem.Tag != dicomtag.PixelData {
			continue
		}

		data := elem.Value[0].(dicom.PixelDataInfo)

		for _, frame := range data.Frames {

			// Encapsulated

			if frame.IsEncapsulated {
				output = append(output, frame.EncapsulatedData.Data)
				continue
			}

			// Unencapsulated

			img := image.NewGray16(image.Rect(0, 0, frame.NativeData.Cols, frame.NativeData.Rows))
			for j := 0; j < len(frame.NativeData.Data); j++ {
				// for now, assume we're not overflowing uint16, assume gray image
				img.SetGray16(j%frame.NativeData.Cols, j/frame.NativeData.Rows, color.Gray16{Y: uint16(frame.NativeData.Data[j][0])})
			}
			buf := new(bytes.Buffer)
			jpeg.Encode(buf, img, &jpeg.Options{Quality: 100})
			output = append(output, buf.Bytes())
		}
	}

	return output, nil
}
