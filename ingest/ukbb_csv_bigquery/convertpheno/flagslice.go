package main

import "strings"

type flagSlice []string

func (i *flagSlice) String() string {
	if i == nil {
		return ""
	}

	return strings.Join([]string(*i), "\t")
}

func (i *flagSlice) Set(value string) error {
	*i = append(*i, value)
	return nil
}
