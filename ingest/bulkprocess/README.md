# Downloading UKBB bulk data (cardiac MRI)
In brief, see [this instruction page](http://biobank.ndph.ox.ac.uk/showcase/instruct/bulk.html). However, these instructions are out of date (give the wrong argument order, etc). See updated instructions below.

Prepare permissions
1. Make sure that you have created a `.ukbkey` file containing the application ID on line 1 and the private key on line 2 (directly downloadable as an attachment from the email that you received from the UKBB). This file should not be readable by anyone without proper UKBB permissions, so consider setting this to be user-readable only.

Download data
1. Download the encrypted file (`ukb21481.enc`) and decrypt it to the encoded file (`ukb21481.enc_ukb`)
1. Extract the list of all samples with the field of interest. 20208 is Heart MRI Long Axis `ukbconv ukb21481.enc_ukb bulk -s20208`
    * Atttempt #2: Try to get all MRI fields at once. `ukbconv ukb21481.enc_ukb bulk -ifields.list`
1. Inspect: `wc -l ukb21481.bulk` and you can see that there is one entry per person for whom this data exists
1. You cannot download more than 1,000 samples' bulk files at a time. So, iteratively do it:
    * For now, just take 50 
    * `head -n 50 ukb21481.bulk > heart.50`
    * `ukbfetch -bheart.50` *(Note: no space between `-b` and `heart.50`)*