# WDL how-to's


## GvsAssignIds.wdl



## GvsImportGenomes.wdl
Currently using the workspace [AoU\_DRC\_WGS\_ah\_ingest\_40k](https://app.terra.bio/#workspaces/allofus-drc-wgs-dev/AoU_DRC_WGS_ah_ingest_40k)

The prerequisites to this step are:

- Reblocking the gvcf for each sampmle
- Assigning a unique gvs id to each sample

We typically want to import genomes into Gvs one sample_set at a time, using the existing sample sets. 

Sometimes the reblocking will fail for a small number of samples, so not all of the samples in the sample set with have the attribute "reblocked\_gvcf\_path". This wdl works on a sample set which means it is passed arrays of attributes all of which need to be filled in. Therefore, we need to create a sample set for loading with the subset of samples that have reblocked gvcfs. To do this, use the notebook `CreateSampleSetForGVSImport-paginated.ipynb`

### CreateSampleSetForGVSImport-paginated.ipynb

1. Get list of samples in sample set

	Currently we are not able to get the list of samples for a sampmle set from the fapi. Instead, go to the data tab, select the sample set you want to process and download the tsv. Unzip the sample set, remove the first line of the membership file (header info) and then extract the sample name column: `cut -f2 sample_set_membership.tsv > samples.txt`


	In the notebook tab, select this notebook and then edit. After the kernel starts, click the jupyter icon to open the filesystem. Navigate to the notebook->edit directory and upload the samples.txt file.
	
2. Edit and run the notebook
	
	Go back to the notebook and set the values in the 3rd box. Set `samples_list_file` to the name of the file you just uploaded. Set `sample_set_name` to the name of the existing sample set you are processing. 
	
	Run each cell of the notebook. Confirm a new sample set with `_gvsload_1` was created.
		
### Run the wdl

Currently the wdl is being imported from the methods repo [here](https://portal.firecloud.org/#methods/andrea_methods/ImportGenomes/68).
