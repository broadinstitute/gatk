<a name="0"></a>
## 0 - Introduction
This page explains the discrepancies between the different "hg19" references.  

There are 4 common "hg19" references, and they are NOT directly interchangeable:
* **hg19** (`ucsc.hg19.fasta`, MD5sum: `a244d8a32473650b25c6e8e1654387d6`)
* **b37** (`Homo_sapiens_assembly19.fasta`, MD5sum: `886ba1559393f75872c1cf459eb57f2d`)
* **GRCh37** (`GRCh37.p13.genome.fasta`, MD5sum: `c140882eb2ea89bc2edfe934d51b66cc`)
* **humanG1Kv37** (`human_g1k_v37.fasta`, MD5sum: `0ce84c872fc0072a885926823dcd0338`)

<a name="0.1"></a>
### 0.1 - Table of Contents
0. [0.0 Introduction](#0)
    1. [0.1 Table of Contents](#0.1)
1. [hg19](#hg19)
2. [GRCh37](#grch37) 
3. [b37](#b37)
4. [HumanG1Kv37](#humanG1Kv37) 
5. [Reference Comparison Table](#comparison)
6. [Additional Information](#additional)

<a name="grch37"></a>
## GRCh37
The Genome Reference Consortium Human Build 37, *GRCh37*, (`GRCh37.p13.genome.fasta`, MD5sum: `c140882eb2ea89bc2edfe934d51b66cc`) is a _Homo Sapiens_ genome reference file built by the Genome Reference Consortium.  This is a baseline human genome reference and serves as the basis for the other three references in this comparison.

For more information on GRCh37, visit the [official Genome Reference Consortium website](https://www.ncbi.nlm.nih.gov/grc).

<a name="hg19"></a>
## hg19
The University of California at Santa Cruz (UCSC) has created a reference based on [GRCh37](#grch37).  This reference is often referred to as *hg19* (`ucsc.hg19.fasta`, MD5sum: `a244d8a32473650b25c6e8e1654387d6`).

This reference contains some alterations from the baseline reference from the Genome Reference Consortium.  These alterations largely consist of contig name changes, however there are *known sequence differences on some contigs as well*.

For details see the [comparison table](#comparison).

<a name="b37"></a>
## b37
The Broad Institute created a human genome reference file based on [GRCh37](#grch37).  This reference is often referred to as *b37* (`Homo_sapiens_assembly19.fasta`, MD5sum: `886ba1559393f75872c1cf459eb57f2d`).

**When people at The Broad Institute's Genomics Platform refer to the _hg19_ reference, they are actually referring to `b37`.**

This reference contains some alterations from the baseline reference from the Genome Reference Consortium.  These alterations largely consist of contig name changes, however there are *known sequence differences on some contigs as well*.  

Anecdotally the changes are for bases for which there was low confidence.  Those low confidence bases were then masked out in the `b37` reference to be the `IUPAC` symbol for any base.  However, _there does not seem to be a detailed comparison readily available_.

For details see the [comparison table](#comparison).

<a name="humanG1Kv37"></a>
## HumanG1Kv37
The *humanG1Kv37* (`human_g1k_v37.fasta`, MD5sum: `0ce84c872fc0072a885926823dcd0338`)  reference is equivalent to [b37](#b37), with the exception that it does not contain the decoy sequence for human herpesvirus 4 type 1 (named _NC_007605_).

<a name="comparison"></a>
## Reference Comparison Table

The specific differences between these four references are detailed in the following table. 

The contigs with identical MD5sums are specified in each row.  In the case that the MD5sum does not match between the references (indicating a sequence difference), the row will have a blank entry for that contig (`----`).

Primary contigs with differing MD5sums are highlighted in red.
Alternate contigs with differing MD5sums are highlighted in orange.

<table>
<tr><th>MD5</th><th>HumanG1Kv37 Contig</th><th>B37 Contig</th><th>HG19 Contig</th><th>GRCh37 Contig</th></tr>
<tr><td>06cbf126247d89664a4faebad130fe9c</td><td>GL000202.1</td><td>GL000202.1</td><td>chr11_gl000202_random</td><td>GL000202.1</td></tr>
<tr><td>0996b4475f353ca98bacb756ac479140</td><td>GL000244.1</td><td>GL000244.1</td><td>chrUn_gl000244</td><td>GL000244.1</td></tr>
<tr><td>118a25ca210cfbcdfb6c2ebb249f9680</td><td>GL000235.1</td><td>GL000235.1</td><td>chrUn_gl000235</td><td>GL000235.1</td></tr>
<tr><td>131b1efc3270cc838686b54e7c34b17b</td><td>GL000238.1</td><td>GL000238.1</td><td>chrUn_gl000238</td><td>GL000238.1</td></tr>
<tr><td>1c1b2cd1fccbc0a99b6a447fa24d1504</td><td>GL000226.1</td><td>GL000226.1</td><td>chrUn_gl000226</td><td>GL000226.1</td></tr>
<tr><td>1d708b54644c26c7e01c2dad5426d38c</td><td>GL000218.1</td><td>GL000218.1</td><td>chrUn_gl000218</td><td>GL000218.1</td></tr>
<tr><td>1d78abec37c15fe29a275eb08d5af236</td><td>GL000249.1</td><td>GL000249.1</td><td>chrUn_gl000249</td><td>GL000249.1</td></tr>
<tr><td>2f8694fc47576bc81b5fe9e7de0ba49e</td><td>GL000242.1</td><td>GL000242.1</td><td>chrUn_gl000242</td><td>GL000242.1</td></tr>
<tr><td>3238fb74ea87ae857f9c7508d315babb</td><td>GL000221.1</td><td>GL000221.1</td><td>chrUn_gl000221</td><td>GL000221.1</td></tr>
<tr><td>325ba9e808f669dfeee210fdd7b470ac</td><td>GL000192.1</td><td>GL000192.1</td><td>chr1_gl000192_random</td><td>GL000192.1</td></tr>
<tr><td>399dfa03bf32022ab52a846f7ca35b30</td><td>GL000223.1</td><td>GL000223.1</td><td>chrUn_gl000223</td><td>GL000223.1</td></tr>
<tr><td>3e06b6741061ad93a8587531307057d8</td><td>GL000232.1</td><td>GL000232.1</td><td>chrUn_gl000232</td><td>GL000232.1</td></tr>
<tr><td>43f69e423533e948bfae5ce1d45bd3f1</td><td>GL000206.1</td><td>GL000206.1</td><td>chr17_gl000206_random</td><td>GL000206.1</td></tr>
<tr><td>445a86173da9f237d7bcf41c6cb8cc62</td><td>GL000240.1</td><td>GL000240.1</td><td>chrUn_gl000240</td><td>GL000240.1</td></tr>
<tr><td>46c2032c37f2ed899eb41c0473319a69</td><td>GL000214.1</td><td>GL000214.1</td><td>chrUn_gl000214</td><td>GL000214.1</td></tr>
<tr><td>563531689f3dbd691331fd6c5730a88b</td><td>GL000212.1</td><td>GL000212.1</td><td>chrUn_gl000212</td><td>GL000212.1</td></tr>
<tr><td>569af3b73522fab4b40995ae4944e78e</td><td>GL000199.1</td><td>GL000199.1</td><td>chr9_gl000199_random</td><td>GL000199.1</td></tr>
<tr><td>5a8e43bec9be36c7b49c84d585107776</td><td>GL000248.1</td><td>GL000248.1</td><td>chrUn_gl000248</td><td>GL000248.1</td></tr>
<tr><td>5d9ec007868d517e73543b005ba48535</td><td>GL000195.1</td><td>GL000195.1</td><td>chr7_gl000195_random</td><td>GL000195.1</td></tr>
<tr><td>5eb3b418480ae67a997957c909375a73</td><td>GL000215.1</td><td>GL000215.1</td><td>chrUn_gl000215</td><td>GL000215.1</td></tr>
<tr><td>63945c3e6962f28ffd469719a747e73c</td><td>GL000225.1</td><td>GL000225.1</td><td>chrUn_gl000225</td><td>GL000225.1</td></tr>
<tr><td>642a232d91c486ac339263820aef7fe0</td><td>GL000216.1</td><td>GL000216.1</td><td>chrUn_gl000216</td><td>GL000216.1</td></tr>
<tr><td>6ac8f815bf8e845bb3031b73f812c012</td><td>GL000194.1</td><td>GL000194.1</td><td>chr4_gl000194_random</td><td>GL000194.1</td></tr>
<tr><td>6d243e18dea1945fb7f2517615b8f52e</td><td>GL000217.1</td><td>GL000217.1</td><td>chrUn_gl000217</td><td>GL000217.1</td></tr>
<tr><td>6f5efdd36643a9b8c8ccad6f2f1edc7b</td><td>GL000197.1</td><td>GL000197.1</td><td>chr8_gl000197_random</td><td>GL000197.1</td></tr>
<tr><td>6fe9abac455169f50470f5a6b01d0f59</td><td>GL000222.1</td><td>GL000222.1</td><td>chrUn_gl000222</td><td>GL000222.1</td></tr>
<tr><td>75e4c8d17cd4addf3917d1703cacaf25</td><td>GL000200.1</td><td>GL000200.1</td><td>chr9_gl000200_random</td><td>GL000200.1</td></tr>
<tr><td>7daaa45c66b288847b9b32b964e623d3</td><td>GL000211.1</td><td>GL000211.1</td><td>chrUn_gl000211</td><td>GL000211.1</td></tr>
<tr><td>7de00226bb7df1c57276ca6baabafd15</td><td>GL000247.1</td><td>GL000247.1</td><td>chrUn_gl000247</td><td>GL000247.1</td></tr>
<tr><td>7e0e2e580297b7764e31dbc80c2540dd</td><td>X</td><td>X</td><td>chrX</td><td>chrX</td></tr>
<tr><td>7fed60298a8d62ff808b74b6ce820001</td><td>GL000233.1</td><td>GL000233.1</td><td>chrUn_gl000233</td><td>GL000233.1</td></tr>
<tr><td>851106a74238044126131ce2a8e5847c</td><td>GL000210.1</td><td>GL000210.1</td><td>chr21_gl000210_random</td><td>GL000210.1</td></tr>
<tr><td>868e7784040da90d900d2d1b667a1383</td><td>GL000198.1</td><td>GL000198.1</td><td>chr9_gl000198_random</td><td>GL000198.1</td></tr>
<tr><td>89bc61960f37d94abf0df2d481ada0ec</td><td>GL000245.1</td><td>GL000245.1</td><td>chrUn_gl000245</td><td>GL000245.1</td></tr>
<tr><td>93f998536b61a56fd0ff47322a911d4b</td><td>GL000234.1</td><td>GL000234.1</td><td>chrUn_gl000234</td><td>GL000234.1</td></tr>
<tr><td>96358c325fe0e70bee73436e8bb14dbd</td><td>GL000203.1</td><td>GL000203.1</td><td>chr17_gl000203_random</td><td>GL000203.1</td></tr>
<tr><td>99795f15702caec4fa1c4e15f8a29c07</td><td>GL000239.1</td><td>GL000239.1</td><td>chrUn_gl000239</td><td>GL000239.1</td></tr>
<tr><td>9d424fdcc98866650b58f004080a992a</td><td>GL000213.1</td><td>GL000213.1</td><td>chrUn_gl000213</td><td>GL000213.1</td></tr>
<tr><td>a4aead23f8053f2655e468bcc6ecdceb</td><td>GL000227.1</td><td>GL000227.1</td><td>chrUn_gl000227</td><td>GL000227.1</td></tr>
<tr><td>aa81be49bf3fe63a79bdc6a6f279abf6</td><td>GL000208.1</td><td>GL000208.1</td><td>chr19_gl000208_random</td><td>GL000208.1</td></tr>
<tr><td>b4eb71ee878d3706246b7c1dbef69299</td><td>GL000230.1</td><td>GL000230.1</td><td>chrUn_gl000230</td><td>GL000230.1</td></tr>
<tr><td>ba8882ce3a1efa2080e5d29b956568a4</td><td>GL000231.1</td><td>GL000231.1</td><td>chrUn_gl000231</td><td>GL000231.1</td></tr>
<tr><td>c5a17c97e2c1a0b6a9cc5a6b064b714f</td><td>GL000228.1</td><td>GL000228.1</td><td>chrUn_gl000228</td><td>GL000228.1</td></tr>
<tr><td>cc34279a7e353136741c9fce79bc4396</td><td>GL000243.1</td><td>GL000243.1</td><td>chrUn_gl000243</td><td>GL000243.1</td></tr>
<tr><td>d0f40ec87de311d8e715b52e4c7062e1</td><td>GL000229.1</td><td>GL000229.1</td><td>chrUn_gl000229</td><td>GL000229.1</td></tr>
<tr><td>d22441398d99caf673e9afb9a1908ec5</td><td>GL000205.1</td><td>GL000205.1</td><td>chr17_gl000205_random</td><td>GL000205.1</td></tr>
<tr><td>d5b2fc04f6b41b212a4198a07f450e20</td><td>GL000224.1</td><td>GL000224.1</td><td>chrUn_gl000224</td><td>GL000224.1</td></tr>
<tr><td>d75b436f50a8214ee9c2a51d30b2c2cc</td><td>GL000191.1</td><td>GL000191.1</td><td>chr1_gl000191_random</td><td>GL000191.1</td></tr>
<tr><td>d92206d1bb4c3b4019c43c0875c06dc0</td><td>GL000196.1</td><td>GL000196.1</td><td>chr8_gl000196_random</td><td>GL000196.1</td></tr>
<tr><td>dbb6e8ece0b5de29da56601613007c2a</td><td>GL000193.1</td><td>GL000193.1</td><td>chr4_gl000193_random</td><td>GL000193.1</td></tr>
<tr><td>dfb7e7ec60ffdcb85cb359ea28454ee9</td><td>GL000201.1</td><td>GL000201.1</td><td>chr9_gl000201_random</td><td>GL000201.1</td></tr>
<tr><td>e0c82e7751df73f4f6d0ed30cdc853c0</td><td>GL000237.1</td><td>GL000237.1</td><td>chrUn_gl000237</td><td>GL000237.1</td></tr>
<tr><td>e4afcd31912af9d9c2546acf1cb23af2</td><td>GL000246.1</td><td>GL000246.1</td><td>chrUn_gl000246</td><td>GL000246.1</td></tr>
<tr><td>ef4258cdc5a45c206cea8fc3e1d858cf</td><td>GL000241.1</td><td>GL000241.1</td><td>chrUn_gl000241</td><td>GL000241.1</td></tr>
<tr><td>efc49c871536fa8d79cb0a06fa739722</td><td>GL000204.1</td><td>GL000204.1</td><td>chr17_gl000204_random</td><td>GL000204.1</td></tr>
<tr><td>f3814841f1939d3ca19072d9e89f3fd7</td><td>GL000207.1</td><td>GL000207.1</td><td>chr18_gl000207_random</td><td>GL000207.1</td></tr>
<tr><td>f40598e2a5a6b26e84a3775e0d1e2c81</td><td>GL000209.1</td><td>GL000209.1</td><td>chr19_gl000209_random</td><td>GL000209.1</td></tr>
<tr><td>f977edd13bac459cb2ed4a5457dba1b3</td><td>GL000219.1</td><td>GL000219.1</td><td>chrUn_gl000219</td><td>GL000219.1</td></tr>
<tr><td>fc35de963c57bf7648429e6454f1c9db</td><td>GL000220.1</td><td>GL000220.1</td><td>chrUn_gl000220</td><td>GL000220.1</td></tr>
<tr><td>fdcd739913efa1fdc64b6c0cd7016779</td><td>GL000236.1</td><td>GL000236.1</td><td>chrUn_gl000236</td><td>GL000236.1</td></tr>
<tr><td>1b22b98cdeb4a9304cb5d48026a85128</td><td>1</td><td>1</td><td>chr1</td><td>chr1</td></tr>
<tr><td>a0d9851da00400dec1098a9255ac712e</td><td>2</td><td>2</td><td>chr2</td><td>chr2</td></tr>
<tr><td>23dccd106897542ad87d2765d28a19a1</td><td>4</td><td>4</td><td>chr4</td><td>chr4</td></tr>
<tr><td>0740173db9ffd264d728f32784845cd7</td><td>5</td><td>5</td><td>chr5</td><td>chr5</td></tr>
<tr><td>1d3a93a248d92a729ee764823acbbc6b</td><td>6</td><td>6</td><td>chr6</td><td>chr6</td></tr>
<tr><td>618366e953d6aaad97dbe4777c29375e</td><td>7</td><td>7</td><td>chr7</td><td>chr7</td></tr>
<tr><td>96f514a9929e410c6651697bded59aec</td><td>8</td><td>8</td><td>chr8</td><td>chr8</td></tr>
<tr><td>3e273117f15e0a400f01055d9f393768</td><td>9</td><td>9</td><td>chr9</td><td>chr9</td></tr>
<tr><td>988c28e000e84c26d552359af1ea2e1d</td><td>10</td><td>10</td><td>chr10</td><td>chr10</td></tr>
<tr><td>98c59049a2df285c76ffb1c6db8f8b96</td><td>11</td><td>11</td><td>chr11</td><td>chr11</td></tr>
<tr><td>51851ac0e1a115847ad36449b0015864</td><td>12</td><td>12</td><td>chr12</td><td>chr12</td></tr>
<tr><td>283f8d7892baa81b510a015719ca7b0b</td><td>13</td><td>13</td><td>chr13</td><td>chr13</td></tr>
<tr><td>98f3cae32b2a2e9524bc19813927542e</td><td>14</td><td>14</td><td>chr14</td><td>chr14</td></tr>
<tr><td>e5645a794a8238215b2cd77acb95a078</td><td>15</td><td>15</td><td>chr15</td><td>chr15</td></tr>
<tr><td>fc9b1a7b42b97a864f56b348b06095e6</td><td>16</td><td>16</td><td>chr16</td><td>chr16</td></tr>
<tr><td>351f64d4f4f9ddd45b35336ad97aa6de</td><td>17</td><td>17</td><td>chr17</td><td>chr17</td></tr>
<tr><td>b15d4b2d29dde9d3e4f93d1d0f2cbc9c</td><td>18</td><td>18</td><td>chr18</td><td>chr18</td></tr>
<tr><td>1aacd71f30db8e561810913e0b72636d</td><td>19</td><td>19</td><td>chr19</td><td>chr19</td></tr>
<tr><td>0dec9660ec1efaaf33281c0d5ea2560f</td><td>20</td><td>20</td><td>chr20</td><td>chr20</td></tr>
<tr><td>2979a6085bfe28e3ad6f552f361ed74d</td><td>21</td><td>21</td><td>chr21</td><td>chr21</td></tr>
<tr><td>a718acaa6135fdca8357d5bfe94211dd</td><td>22</td><td>22</td><td>chr22</td><td>chr22</td></tr>
<tr bgcolor="#FD1010"><td>1fa3474750af0948bdf97d5a0ee52e51</td><td>Y</td><td>Y</td><td>----</td><td>----</td></tr>
<tr bgcolor="#EB8B31"><td>6743bd63b3ff2b5b8985d8933c53290a</td><td>----</td><td>NC_007605</td><td>----</td><td>----</td></tr>
<tr bgcolor="#FD1010"><td>c68f52674c9fb33aef52dcf399755519</td><td>MT</td><td>MT</td><td>----</td><td>chrM</td></tr>
<tr bgcolor="#FD1010"><td>fdfd811849cc2fadebc929bb925902e5</td><td>3</td><td>3</td><td>----</td><td>----</td></tr>
<tr bgcolor="#EB8B31"><td>094d037050cad692b57ea12c4fef790f</td><td>----</td><td>----</td><td>chr6_qbl_hap6</td><td>GL000255.1</td></tr>
<tr bgcolor="#EB8B31"><td>18c17e1641ef04873b15f40f6c8659a4</td><td>----</td><td>----</td><td>chr6_cox_hap2</td><td>GL000251.1</td></tr>
<tr bgcolor="#FD1010"><td>1e86411d73e6f00a10590f976be01623</td><td>----</td><td>----</td><td>chrY</td><td>chrY</td></tr>
<tr bgcolor="#EB8B31"><td>2a3c677c426a10e137883ae1ffb8da3f</td><td>----</td><td>----</td><td>chr6_dbb_hap3</td><td>GL000252.1</td></tr>
<tr bgcolor="#EB8B31"><td>3b6d666200e72bcc036bf88a4d7e0749</td><td>----</td><td>----</td><td>chr6_ssto_hap7</td><td>GL000256.1</td></tr>
<tr bgcolor="#FD1010"><td>641e4338fa8d52a5b781bd2a2c08d3c3</td><td>----</td><td>----</td><td>chr3</td><td>chr3</td></tr>
<tr bgcolor="#EB8B31"><td>9d51d4152174461cd6715c7ddc588dc8</td><td>----</td><td>----</td><td>chr6_mann_hap4</td><td>GL000253.1</td></tr>
<tr bgcolor="#FD1010"><td>d2ed829b8a1628d16cbeee88e88e39eb</td><td>----</td><td>----</td><td>chrM</td><td>----</td></tr>
<tr bgcolor="#EB8B31"><td>d89517b400226d3b56e753972a7cad67</td><td>----</td><td>----</td><td>chr17_ctg5_hap1</td><td>GL000258.1</td></tr>
<tr bgcolor="#EB8B31"><td>efed415dd8742349cb7aaca054675b9a</td><td>----</td><td>----</td><td>chr6_mcf_hap5</td><td>GL000254.1</td></tr>
<tr bgcolor="#EB8B31"><td>fa24f81b680df26bcfb6d69b784fbe36</td><td>----</td><td>----</td><td>chr4_ctg9_hap1</td><td>GL000257.1</td></tr>
<tr bgcolor="#EB8B31"><td>fe71bc63420d666884f37a3ad79f3317</td><td>----</td><td>----</td><td>chr6_apd_hap1</td><td>GL000250.1</td></tr>
<tr bgcolor="#EB8B31"><td>0386df1d3476e6649f919195cc072fc7</td><td>----</td><td>----</td><td>----</td><td>GL383574.1</td></tr>
<tr bgcolor="#EB8B31"><td>03de7a950b56720768373120bbddf693</td><td>----</td><td>----</td><td>----</td><td>GL383558.1</td></tr>
<tr bgcolor="#EB8B31"><td>03f5fa89e52d0fe155d2e3968bf2eeb7</td><td>----</td><td>----</td><td>----</td><td>GL339450.1</td></tr>
<tr bgcolor="#EB8B31"><td>049056a72b5aee0b3f876ddf554f0208</td><td>----</td><td>----</td><td>----</td><td>JH720447.1</td></tr>
<tr bgcolor="#EB8B31"><td>063358c8e7f81361b959efab7b3f15cc</td><td>----</td><td>----</td><td>----</td><td>GL383565.1</td></tr>
<tr bgcolor="#EB8B31"><td>07f56906bd56829f146dc0bf4b158603</td><td>----</td><td>----</td><td>----</td><td>GL383538.1</td></tr>
<tr bgcolor="#EB8B31"><td>09ce2d45f1f973e347e22ad1e3cf06fb</td><td>----</td><td>----</td><td>----</td><td>GL383568.1</td></tr>
<tr bgcolor="#EB8B31"><td>09d4cb1070e1c521d6e86e7038824c1c</td><td>----</td><td>----</td><td>----</td><td>GL383544.1</td></tr>
<tr bgcolor="#EB8B31"><td>0aee7c3e4bcc4c942230508c7836069b</td><td>----</td><td>----</td><td>----</td><td>GL582967.1</td></tr>
<tr bgcolor="#EB8B31"><td>0bdfd2a40e1ceab32d71d6f1c9a6ca32</td><td>----</td><td>----</td><td>----</td><td>JH636056.1</td></tr>
<tr bgcolor="#EB8B31"><td>0c787911df2449cbba8609bebf897ecb</td><td>----</td><td>----</td><td>----</td><td>GL383541.1</td></tr>
<tr bgcolor="#EB8B31"><td>0d12851232bcd8250e6dd61e3e7fd6a2</td><td>----</td><td>----</td><td>----</td><td>JH636054.1</td></tr>
<tr bgcolor="#EB8B31"><td>0db85f8e0ff66a470b46801c9892e471</td><td>----</td><td>----</td><td>----</td><td>JH806589.1</td></tr>
<tr bgcolor="#EB8B31"><td>0de582e28ae8c127d978b43e12a4f499</td><td>----</td><td>----</td><td>----</td><td>KE332500.1</td></tr>
<tr bgcolor="#EB8B31"><td>0f0364ed52ebe7757feea96ce623239f</td><td>----</td><td>----</td><td>----</td><td>GL383529.1</td></tr>
<tr bgcolor="#EB8B31"><td>10b523cdfd4f3707276ec92f0f9cddfb</td><td>----</td><td>----</td><td>----</td><td>JH159133.1</td></tr>
<tr bgcolor="#EB8B31"><td>1135da7b213739cfb3bdf2741c0c8083</td><td>----</td><td>----</td><td>----</td><td>KE332506.1</td></tr>
<tr bgcolor="#EB8B31"><td>122be4e189778434d8845fd5fd2c9a6b</td><td>----</td><td>----</td><td>----</td><td>GL949748.1</td></tr>
<tr bgcolor="#EB8B31"><td>12406aad3f3da31bda9c21a1aa0e16b6</td><td>----</td><td>----</td><td>----</td><td>GL383539.1</td></tr>
<tr bgcolor="#EB8B31"><td>12a3180640a49f33c960eb12ca61a6c4</td><td>----</td><td>----</td><td>----</td><td>GL383542.1</td></tr>
<tr bgcolor="#EB8B31"><td>13a8dc0d93c1bf1ae397593eba841721</td><td>----</td><td>----</td><td>----</td><td>JH806573.1</td></tr>
<tr bgcolor="#EB8B31"><td>1452be48789c27311d94561610f6d5af</td><td>----</td><td>----</td><td>----</td><td>JH159132.1</td></tr>
<tr bgcolor="#EB8B31"><td>157eecba5817aa1781a7bc4a9b60f933</td><td>----</td><td>----</td><td>----</td><td>JH720446.1</td></tr>
<tr bgcolor="#EB8B31"><td>15a2182cf9d3a55a7809adadc0775e03</td><td>----</td><td>----</td><td>----</td><td>JH636060.1</td></tr>
<tr bgcolor="#EB8B31"><td>16197ace4bfafdc2354857b98fc2a794</td><td>----</td><td>----</td><td>----</td><td>GL582971.1</td></tr>
<tr bgcolor="#EB8B31"><td>1668b0eb03be297f66837b46b5b73ac7</td><td>----</td><td>----</td><td>----</td><td>JH806600.2</td></tr>
<tr bgcolor="#EB8B31"><td>16a9ef53c176dcd2cf029940cbc29382</td><td>----</td><td>----</td><td>----</td><td>JH720452.1</td></tr>
<tr bgcolor="#EB8B31"><td>18c01f6e62136005ce1b2f2f33173f02</td><td>----</td><td>----</td><td>----</td><td>GL582974.1</td></tr>
<tr bgcolor="#EB8B31"><td>18fd9605f12ec0982adcf9e908f53331</td><td>----</td><td>----</td><td>----</td><td>GL383523.1</td></tr>
<tr bgcolor="#EB8B31"><td>190b4d20cab29ddba5f2495a7117f7cf</td><td>----</td><td>----</td><td>----</td><td>JH636055.1</td></tr>
<tr bgcolor="#EB8B31"><td>1919a95f3ea48fde56ba925295086028</td><td>----</td><td>----</td><td>----</td><td>GL383582.2</td></tr>
<tr bgcolor="#EB8B31"><td>192dead6bc331a0dcbd1ba9d3d8a6f80</td><td>----</td><td>----</td><td>----</td><td>JH806581.1</td></tr>
<tr bgcolor="#EB8B31"><td>1b6fa375fdf382778e6645d822d12254</td><td>----</td><td>----</td><td>----</td><td>GL339449.2</td></tr>
<tr bgcolor="#EB8B31"><td>1d217e666c48ef7f2ec81946f4f4bfcc</td><td>----</td><td>----</td><td>----</td><td>GL383521.1</td></tr>
<tr bgcolor="#EB8B31"><td>1d2933992f087a832718c9de19d4ceab</td><td>----</td><td>----</td><td>----</td><td>KE332499.1</td></tr>
<tr bgcolor="#EB8B31"><td>20d5046bbd2a21729fdd64fa94bdd5a1</td><td>----</td><td>----</td><td>----</td><td>GL949742.1</td></tr>
<tr bgcolor="#EB8B31"><td>20da91baf79b2e14b605a8ebe1f3704e</td><td>----</td><td>----</td><td>----</td><td>GL877872.1</td></tr>
<tr bgcolor="#EB8B31"><td>2118ff7bca8f75acc4629ab88bae1c2e</td><td>----</td><td>----</td><td>----</td><td>GL582968.1</td></tr>
<tr bgcolor="#EB8B31"><td>23aea04f46682e2a2be1a5ff3934a9fe</td><td>----</td><td>----</td><td>----</td><td>GL383540.1</td></tr>
<tr bgcolor="#EB8B31"><td>24f6ccbfc261e62451042a9713be6280</td><td>----</td><td>----</td><td>----</td><td>JH636059.1</td></tr>
<tr bgcolor="#EB8B31"><td>2536c2286fdeb98404ac410dbf528a3e</td><td>----</td><td>----</td><td>----</td><td>GL949745.1</td></tr>
<tr bgcolor="#EB8B31"><td>260aca5e1ff29b6ed3d3cd9e438f4219</td><td>----</td><td>----</td><td>----</td><td>JH636053.3</td></tr>
<tr bgcolor="#EB8B31"><td>2948653361f974fbed3e26a4dfbf332c</td><td>----</td><td>----</td><td>----</td><td>GL383528.1</td></tr>
<tr bgcolor="#EB8B31"><td>2cfa9ec8f70be88f95411dac6efb24c1</td><td>----</td><td>----</td><td>----</td><td>JH159141.2</td></tr>
<tr bgcolor="#EB8B31"><td>2e0bec27cfa9b440c746be52187fab0b</td><td>----</td><td>----</td><td>----</td><td>GL383560.1</td></tr>
<tr bgcolor="#EB8B31"><td>2f4f58e3b3a95bed1132833156340778</td><td>----</td><td>----</td><td>----</td><td>GL383555.1</td></tr>
<tr bgcolor="#EB8B31"><td>2fc316247e162f76a01012bbd9b665e6</td><td>----</td><td>----</td><td>----</td><td>JH806603.1</td></tr>
<tr bgcolor="#EB8B31"><td>321b3324431ae40e90a4117ecc07e93c</td><td>----</td><td>----</td><td>----</td><td>JH591185.1</td></tr>
<tr bgcolor="#EB8B31"><td>32ceefe714becfd36f207c5bffca4ba7</td><td>----</td><td>----</td><td>----</td><td>GL877871.1</td></tr>
<tr bgcolor="#EB8B31"><td>3342f6c21fa2dc364925712e0d52ed2d</td><td>----</td><td>----</td><td>----</td><td>JH159148.1</td></tr>
<tr bgcolor="#EB8B31"><td>348782844074360bdc8b6416e16cc5d0</td><td>----</td><td>----</td><td>----</td><td>GL383549.1</td></tr>
<tr bgcolor="#EB8B31"><td>349e96f115f829409bd1087b5fb684ca</td><td>----</td><td>----</td><td>----</td><td>GL383519.1</td></tr>
<tr bgcolor="#EB8B31"><td>35889722e6212fc9499e06e630268101</td><td>----</td><td>----</td><td>----</td><td>JH806588.1</td></tr>
<tr bgcolor="#EB8B31"><td>369f03e72d44461eab4542c58f3b5dcc</td><td>----</td><td>----</td><td>----</td><td>KE332501.1</td></tr>
<tr bgcolor="#EB8B31"><td>384b5b32f0ea2cfd15ac268a2ce07909</td><td>----</td><td>----</td><td>----</td><td>JH159146.1</td></tr>
<tr bgcolor="#EB8B31"><td>3893d44dfad7ce35744b2bde1e43bbd3</td><td>----</td><td>----</td><td>----</td><td>GL383576.1</td></tr>
<tr bgcolor="#EB8B31"><td>38e72cd57edb75d967ac2613d61d297d</td><td>----</td><td>----</td><td>----</td><td>GL383563.2</td></tr>
<tr bgcolor="#EB8B31"><td>3b40b7fdb005a1ce00efaa3310148852</td><td>----</td><td>----</td><td>----</td><td>KE332502.1</td></tr>
<tr bgcolor="#EB8B31"><td>3c5f20fb0744b7658d37d4ed79a286d1</td><td>----</td><td>----</td><td>----</td><td>GL383520.1</td></tr>
<tr bgcolor="#EB8B31"><td>3dd30a7638c3a3c518fc15571546b1be</td><td>----</td><td>----</td><td>----</td><td>GL877875.1</td></tr>
<tr bgcolor="#EB8B31"><td>3e0825dd23c9fce74a88d863e33c42b7</td><td>----</td><td>----</td><td>----</td><td>JH159149.1</td></tr>
<tr bgcolor="#EB8B31"><td>40015159c7da8f06875bb558587e3f07</td><td>----</td><td>----</td><td>----</td><td>GL383571.1</td></tr>
<tr bgcolor="#EB8B31"><td>404580d8ad56ded0fb33642c8b99c28b</td><td>----</td><td>----</td><td>----</td><td>GL383537.1</td></tr>
<tr bgcolor="#EB8B31"><td>4112dee892050e18ad279b8ebdcc5d48</td><td>----</td><td>----</td><td>----</td><td>GL383578.1</td></tr>
<tr bgcolor="#EB8B31"><td>41cf432193561894813a30da6e682e5b</td><td>----</td><td>----</td><td>----</td><td>KB021647.1</td></tr>
<tr bgcolor="#EB8B31"><td>447fe0ff3103170150280c775095eebf</td><td>----</td><td>----</td><td>----</td><td>JH806593.1</td></tr>
<tr bgcolor="#EB8B31"><td>44d5da56e5ec6ae0b9ebd354e9b47cfa</td><td>----</td><td>----</td><td>----</td><td>JH159150.3</td></tr>
<tr bgcolor="#EB8B31"><td>474baa8f6c6684c55bbc2a10bfa84baf</td><td>----</td><td>----</td><td>----</td><td>GL949747.1</td></tr>
<tr bgcolor="#EB8B31"><td>4791ba11d768da2cc1346d37a558047a</td><td>----</td><td>----</td><td>----</td><td>GL582972.1</td></tr>
<tr bgcolor="#EB8B31"><td>485c442c93fe19514153702f0c84d952</td><td>----</td><td>----</td><td>----</td><td>GL582966.2</td></tr>
<tr bgcolor="#EB8B31"><td>4a3d54bda53308ca941d6d0e794b05cb</td><td>----</td><td>----</td><td>----</td><td>GL383554.1</td></tr>
<tr bgcolor="#EB8B31"><td>4ad67c3a4e85f8b2cd54a8ef2aab4426</td><td>----</td><td>----</td><td>----</td><td>GL383557.1</td></tr>
<tr bgcolor="#EB8B31"><td>4bc4f02a4fca2c9d70646455bee8066e</td><td>----</td><td>----</td><td>----</td><td>JH806596.1</td></tr>
<tr bgcolor="#EB8B31"><td>50fd52ddb8ad2b024fb8b83a5c90a642</td><td>----</td><td>----</td><td>----</td><td>JH806577.1</td></tr>
<tr bgcolor="#EB8B31"><td>5404455aab275489bc8e6c9fb3ead5cb</td><td>----</td><td>----</td><td>----</td><td>GL383573.1</td></tr>
<tr bgcolor="#EB8B31"><td>54abe159678a84e88ceb2d5271027628</td><td>----</td><td>----</td><td>----</td><td>KB663609.1</td></tr>
<tr bgcolor="#EB8B31"><td>5835d9de56b65cefb9406d104d64531e</td><td>----</td><td>----</td><td>----</td><td>GL383536.1</td></tr>
<tr bgcolor="#EB8B31"><td>5950c02594cedbdf0fea5e8335e7cf80</td><td>----</td><td>----</td><td>----</td><td>JH591181.2</td></tr>
<tr bgcolor="#EB8B31"><td>5b90c3ac4e5938b400fcc2c29f3017bc</td><td>----</td><td>----</td><td>----</td><td>GL949744.1</td></tr>
<tr bgcolor="#EB8B31"><td>5b9d9fb059071e552bba531b81bd3472</td><td>----</td><td>----</td><td>----</td><td>KB663606.1</td></tr>
<tr bgcolor="#EB8B31"><td>5c3a364520bf7ed46894abdce8f6e032</td><td>----</td><td>----</td><td>----</td><td>GL877876.1</td></tr>
<tr bgcolor="#EB8B31"><td>5eb6da458990f121fae13ff83a4bcbca</td><td>----</td><td>----</td><td>----</td><td>JH806576.1</td></tr>
<tr bgcolor="#EB8B31"><td>5f9dc3f86463d08a1383cca5f285b7ad</td><td>----</td><td>----</td><td>----</td><td>GL383570.1</td></tr>
<tr bgcolor="#EB8B31"><td>5fae03628eb9a445571bac107823b394</td><td>----</td><td>----</td><td>----</td><td>KE332505.1</td></tr>
<tr bgcolor="#EB8B31"><td>620913159e2fbd4e931ac120e3c584c9</td><td>----</td><td>----</td><td>----</td><td>GL383526.1</td></tr>
<tr bgcolor="#EB8B31"><td>659b65783878ace88f4c4b165f239363</td><td>----</td><td>----</td><td>----</td><td>JH636052.4</td></tr>
<tr bgcolor="#EB8B31"><td>675046e52613269a7c2e803525bb5a33</td><td>----</td><td>----</td><td>----</td><td>JH720445.1</td></tr>
<tr bgcolor="#EB8B31"><td>67f26a755ca4c6ca9a8f567d80d15fb9</td><td>----</td><td>----</td><td>----</td><td>JH159131.1</td></tr>
<tr bgcolor="#EB8B31"><td>68391fb8f16a37b63f607b76702de3b1</td><td>----</td><td>----</td><td>----</td><td>GL383583.1</td></tr>
<tr bgcolor="#EB8B31"><td>69490aa24b00717f2b11c095a5339516</td><td>----</td><td>----</td><td>----</td><td>JH720454.3</td></tr>
<tr bgcolor="#EB8B31"><td>6b604cf3e324680b72716e814d805944</td><td>----</td><td>----</td><td>----</td><td>GL383530.1</td></tr>
<tr bgcolor="#EB8B31"><td>6b862a953dfe724a1f48eaf12a3b948a</td><td>----</td><td>----</td><td>----</td><td>JH806578.1</td></tr>
<tr bgcolor="#EB8B31"><td>6c22616c927261b8e5fc90028c780f00</td><td>----</td><td>----</td><td>----</td><td>JH720444.2</td></tr>
<tr bgcolor="#EB8B31"><td>6cba57c0e509ab785d3869134979b668</td><td>----</td><td>----</td><td>----</td><td>GL383548.1</td></tr>
<tr bgcolor="#EB8B31"><td>6d728406957c5c7fb158dbdb7efef2b7</td><td>----</td><td>----</td><td>----</td><td>GL383527.1</td></tr>
<tr bgcolor="#EB8B31"><td>6d85d704338ba29941aca4d278c7eb4a</td><td>----</td><td>----</td><td>----</td><td>KB663608.1</td></tr>
<tr bgcolor="#EB8B31"><td>6fbc7007a5ff8aae8b28ca52dd6d5571</td><td>----</td><td>----</td><td>----</td><td>GL383531.1</td></tr>
<tr bgcolor="#EB8B31"><td>71a0d13c09c3e7ee64c7740e1425f20e</td><td>----</td><td>----</td><td>----</td><td>JH591183.1</td></tr>
<tr bgcolor="#EB8B31"><td>73b240dd73b8bddcab281e265c9d759a</td><td>----</td><td>----</td><td>----</td><td>GL383559.2</td></tr>
<tr bgcolor="#EB8B31"><td>73d39b5d51e6e2e8d9549bb85d7dae04</td><td>----</td><td>----</td><td>----</td><td>JH806579.1</td></tr>
<tr bgcolor="#EB8B31"><td>741179e4ee12c60fbcc6eba4a5c7695b</td><td>----</td><td>----</td><td>----</td><td>JH806598.1</td></tr>
<tr bgcolor="#EB8B31"><td>74cff045a9cd92b7f571a756f248d16a</td><td>----</td><td>----</td><td>----</td><td>GL383524.1</td></tr>
<tr bgcolor="#EB8B31"><td>7b556f03729e304a286c8d7ef0f0c10e</td><td>----</td><td>----</td><td>----</td><td>GL383547.1</td></tr>
<tr bgcolor="#EB8B31"><td>7b6d6d01c18e91fc07f727ade2450f46</td><td>----</td><td>----</td><td>----</td><td>JH806580.1</td></tr>
<tr bgcolor="#EB8B31"><td>7d007a35ff02e56325881c68bb17b565</td><td>----</td><td>----</td><td>----</td><td>GL949752.1</td></tr>
<tr bgcolor="#EB8B31"><td>7e0afbdc97540aa0b101228b7bd331fb</td><td>----</td><td>----</td><td>----</td><td>JH806591.1</td></tr>
<tr bgcolor="#EB8B31"><td>8213c58e2c1c22397f0ad9d0d901bbdf</td><td>----</td><td>----</td><td>----</td><td>JH159135.2</td></tr>
<tr bgcolor="#EB8B31"><td>856a46516332f58a35eeb4f84d17febc</td><td>----</td><td>----</td><td>----</td><td>JH159142.2</td></tr>
<tr bgcolor="#EB8B31"><td>883b29a1e5975e0f3139c183fbe2596d</td><td>----</td><td>----</td><td>----</td><td>GL383581.1</td></tr>
<tr bgcolor="#EB8B31"><td>8a92722deabdf885d1aebfa8881d5903</td><td>----</td><td>----</td><td>----</td><td>GL949753.1</td></tr>
<tr bgcolor="#EB8B31"><td>8ac2dc8046e4bd0d6d46e827ff05ecd1</td><td>----</td><td>----</td><td>----</td><td>JH806575.1</td></tr>
<tr bgcolor="#EB8B31"><td>8ac9fb9d942dba38bfd30f8d767f4bba</td><td>----</td><td>----</td><td>----</td><td>JH159136.1</td></tr>
<tr bgcolor="#EB8B31"><td>8b1d46e46d3083625eac92e9363773dd</td><td>----</td><td>----</td><td>----</td><td>KB663607.2</td></tr>
<tr bgcolor="#EB8B31"><td>8d13c3e7cbb2b7e1a3225c5a54fe8f44</td><td>----</td><td>----</td><td>----</td><td>GL383569.1</td></tr>
<tr bgcolor="#EB8B31"><td>8e1004755b0574b2f855130c943fbd8e</td><td>----</td><td>----</td><td>----</td><td>KB663603.1</td></tr>
<tr bgcolor="#EB8B31"><td>8e4f862d5b37504199902c6685b7fee5</td><td>----</td><td>----</td><td>----</td><td>GL383533.1</td></tr>
<tr bgcolor="#EB8B31"><td>8ede6eec21d781c22a9801f51433fcd6</td><td>----</td><td>----</td><td>----</td><td>JH806584.1</td></tr>
<tr bgcolor="#EB8B31"><td>8fc7aaa775b43df3d77c9782a140a981</td><td>----</td><td>----</td><td>----</td><td>GL383572.1</td></tr>
<tr bgcolor="#EB8B31"><td>902d62224f09e59cb9c6c44f71b5fca3</td><td>----</td><td>----</td><td>----</td><td>GL383580.1</td></tr>
<tr bgcolor="#EB8B31"><td>90ad438579d919fd20c42bb4f48de64b</td><td>----</td><td>----</td><td>----</td><td>JH591184.1</td></tr>
<tr bgcolor="#EB8B31"><td>9133580f75d0ffa745af12953d65a4db</td><td>----</td><td>----</td><td>----</td><td>GL383550.1</td></tr>
<tr bgcolor="#EB8B31"><td>93001afcfc8594885490513c4ffe243e</td><td>----</td><td>----</td><td>----</td><td>GL582977.2</td></tr>
<tr bgcolor="#EB8B31"><td>93a798f03267e553445c7456c6f7ee49</td><td>----</td><td>----</td><td>----</td><td>JH159143.1</td></tr>
<tr bgcolor="#EB8B31"><td>94409f94ca59e67f811cd36ab133a82c</td><td>----</td><td>----</td><td>----</td><td>GL877873.1</td></tr>
<tr bgcolor="#EB8B31"><td>955e16dfcdb2d28a334349dfa39f2ed4</td><td>----</td><td>----</td><td>----</td><td>KB021645.1</td></tr>
<tr bgcolor="#EB8B31"><td>976c3a7c4051dd9ce879833f4a764289</td><td>----</td><td>----</td><td>----</td><td>JH806590.2</td></tr>
<tr bgcolor="#EB8B31"><td>978987018f1a910273ebcc387e038de8</td><td>----</td><td>----</td><td>----</td><td>GL383518.1</td></tr>
<tr bgcolor="#EB8B31"><td>9bb3fbcd1fc9c35884e0987755c55667</td><td>----</td><td>----</td><td>----</td><td>GL949741.1</td></tr>
<tr bgcolor="#EB8B31"><td>9bed9883963242dd74b883218d5f17bf</td><td>----</td><td>----</td><td>----</td><td>GL383532.1</td></tr>
<tr bgcolor="#EB8B31"><td>9d197695e8a47d4c30c891a53a0fd588</td><td>----</td><td>----</td><td>----</td><td>JH720451.1</td></tr>
<tr bgcolor="#EB8B31"><td>9dcedb7219aa23057244ca9a446f01ac</td><td>----</td><td>----</td><td>----</td><td>JH806592.1</td></tr>
<tr bgcolor="#EB8B31"><td>9e1fc7ed55646756ce109b12b82ff192</td><td>----</td><td>----</td><td>----</td><td>JH159147.1</td></tr>
<tr bgcolor="#EB8B31"><td>a009cf3116a844d7b2d467e672931bc5</td><td>----</td><td>----</td><td>----</td><td>GL383564.1</td></tr>
<tr bgcolor="#EB8B31"><td>a0584071d5a8e88fda38d4cca38704cb</td><td>----</td><td>----</td><td>----</td><td>KB663605.1</td></tr>
<tr bgcolor="#EB8B31"><td>a0bce2b33eb96adcb750622527225e7d</td><td>----</td><td>----</td><td>----</td><td>JH720443.2</td></tr>
<tr bgcolor="#EB8B31"><td>a0f25165c6537c9861cc1231f710e99f</td><td>----</td><td>----</td><td>----</td><td>GL383566.1</td></tr>
<tr bgcolor="#EB8B31"><td>a1dab5e9bbedd3539ace29af1f9d6139</td><td>----</td><td>----</td><td>----</td><td>JH806595.1</td></tr>
<tr bgcolor="#EB8B31"><td>a260ca7327d292deefef4f5fc7346dc4</td><td>----</td><td>----</td><td>----</td><td>GL383561.2</td></tr>
<tr bgcolor="#EB8B31"><td>a2ecd2eb53eb1737423d5a637e4374a9</td><td>----</td><td>----</td><td>----</td><td>JH636058.1</td></tr>
<tr bgcolor="#EB8B31"><td>a3bf927c2422ea0a661640669efd1081</td><td>----</td><td>----</td><td>----</td><td>GL383577.1</td></tr>
<tr bgcolor="#EB8B31"><td>a4053747fc0cf1e03fa6ae9cd5f821d0</td><td>----</td><td>----</td><td>----</td><td>KB021648.1</td></tr>
<tr bgcolor="#EB8B31"><td>aa5b0a15acec3c6177db764bd103d8a0</td><td>----</td><td>----</td><td>----</td><td>JH720453.1</td></tr>
<tr bgcolor="#EB8B31"><td>ab73a8d586ef4fc44dd063730b6aef39</td><td>----</td><td>----</td><td>----</td><td>GL582973.1</td></tr>
<tr bgcolor="#EB8B31"><td>abb9297c8b9dfc3013d416c803ff486c</td><td>----</td><td>----</td><td>----</td><td>JH806601.1</td></tr>
<tr bgcolor="#EB8B31"><td>ac9c384b2fc322b684128f1baf75785e</td><td>----</td><td>----</td><td>----</td><td>JH591182.1</td></tr>
<tr bgcolor="#EB8B31"><td>adb23c033121d433739de02cfa00c9fb</td><td>----</td><td>----</td><td>----</td><td>JH806582.2</td></tr>
<tr bgcolor="#EB8B31"><td>adec63ae44a39d716808cfee03b7a870</td><td>----</td><td>----</td><td>----</td><td>JH806599.1</td></tr>
<tr bgcolor="#EB8B31"><td>afb0d13ed9fa7518989caa0ec55aeb96</td><td>----</td><td>----</td><td>----</td><td>GL949750.1</td></tr>
<tr bgcolor="#EB8B31"><td>b293c854ddcbc316cb1d449bca46fbb3</td><td>----</td><td>----</td><td>----</td><td>JH159137.1</td></tr>
<tr bgcolor="#EB8B31"><td>b8864877618b25fc14f80e8538f23b77</td><td>----</td><td>----</td><td>----</td><td>GL949751.1</td></tr>
<tr bgcolor="#EB8B31"><td>b96f5e6bc844e8392d4e442aa7557e15</td><td>----</td><td>----</td><td>----</td><td>KE332495.1</td></tr>
<tr bgcolor="#EB8B31"><td>ba6a3b1599661e674918200a8d1333d3</td><td>----</td><td>----</td><td>----</td><td>JH806594.1</td></tr>
<tr bgcolor="#EB8B31"><td>bc6f64b0c4c934c2cea52bbe98639c79</td><td>----</td><td>----</td><td>----</td><td>JH806583.1</td></tr>
<tr bgcolor="#EB8B31"><td>bc79d1abee7076ea672293e12bd7ccb9</td><td>----</td><td>----</td><td>----</td><td>GL383517.1</td></tr>
<tr bgcolor="#EB8B31"><td>bd742a610e4bbc28fc00aaf71dfdc15d</td><td>----</td><td>----</td><td>----</td><td>KB021644.1</td></tr>
<tr bgcolor="#EB8B31"><td>be51fd8c00d62c3efc077a8e882062a4</td><td>----</td><td>----</td><td>----</td><td>JH806586.1</td></tr>
<tr bgcolor="#EB8B31"><td>bed6a2667e8452a176e93e921e0c21f6</td><td>----</td><td>----</td><td>----</td><td>GL383556.1</td></tr>
<tr bgcolor="#EB8B31"><td>c27dc6fea378fecf178a44682257c25e</td><td>----</td><td>----</td><td>----</td><td>GL383545.1</td></tr>
<tr bgcolor="#EB8B31"><td>c28f12c6ee0dec4cc6995766a710960c</td><td>----</td><td>----</td><td>----</td><td>GL383552.1</td></tr>
<tr bgcolor="#EB8B31"><td>c6ff49147dedce02366d6ade10580611</td><td>----</td><td>----</td><td>----</td><td>GL383534.2</td></tr>
<tr bgcolor="#EB8B31"><td>c86ffa095c924372aa455e43e61c96e8</td><td>----</td><td>----</td><td>----</td><td>JH159145.1</td></tr>
<tr bgcolor="#EB8B31"><td>ca0e3270f27bbee944844e44ec76659d</td><td>----</td><td>----</td><td>----</td><td>JH806602.1</td></tr>
<tr bgcolor="#EB8B31"><td>caebc01e3f44f7b2a559179b0261b77e</td><td>----</td><td>----</td><td>----</td><td>GL383535.1</td></tr>
<tr bgcolor="#EB8B31"><td>cca1c60136ec678eeef374134cd07a90</td><td>----</td><td>----</td><td>----</td><td>GL877877.2</td></tr>
<tr bgcolor="#EB8B31"><td>cdab95f32513753b3c0add3014afad3b</td><td>----</td><td>----</td><td>----</td><td>JH806587.1</td></tr>
<tr bgcolor="#EB8B31"><td>d08cc284ad35f0bd1eafb443c23ad8bd</td><td>----</td><td>----</td><td>----</td><td>JH159138.1</td></tr>
<tr bgcolor="#EB8B31"><td>d0b63f9cef6c4d382e49636465eab851</td><td>----</td><td>----</td><td>----</td><td>GL877870.2</td></tr>
<tr bgcolor="#EB8B31"><td>d0caa7bf982cf1e6ca8c8b833f56a21c</td><td>----</td><td>----</td><td>----</td><td>JH806597.1</td></tr>
<tr bgcolor="#EB8B31"><td>d4e2cf05984db16a78c953b898f5a86e</td><td>----</td><td>----</td><td>----</td><td>GL582969.1</td></tr>
<tr bgcolor="#EB8B31"><td>d76e635e75bc038782fb3d0c195d33fb</td><td>----</td><td>----</td><td>----</td><td>GL949746.1</td></tr>
<tr bgcolor="#EB8B31"><td>d8ef242a7373ff5657c8311b92dabfde</td><td>----</td><td>----</td><td>----</td><td>JH591186.1</td></tr>
<tr bgcolor="#EB8B31"><td>d9015dd9a0916a98ed8ab99fd3cdd012</td><td>----</td><td>----</td><td>----</td><td>GL383567.1</td></tr>
<tr bgcolor="#EB8B31"><td>d96719c32333013a51c4d6d3261f984f</td><td>----</td><td>----</td><td>----</td><td>GL383551.1</td></tr>
<tr bgcolor="#EB8B31"><td>d97cf75e24ed1370388fedf523faa7ab</td><td>----</td><td>----</td><td>----</td><td>KE332496.1</td></tr>
<tr bgcolor="#EB8B31"><td>da648c938f1bb43b41d254bd9a015cfb</td><td>----</td><td>----</td><td>----</td><td>JH159139.1</td></tr>
<tr bgcolor="#EB8B31"><td>dada6dd12ec844a3a13f547f4946428e</td><td>----</td><td>----</td><td>----</td><td>JH159140.1</td></tr>
<tr bgcolor="#EB8B31"><td>dd0bc538e31f35af2073daec1f378147</td><td>----</td><td>----</td><td>----</td><td>JH159144.1</td></tr>
<tr bgcolor="#EB8B31"><td>dd784bb8074d6f5b949464ffea8c6901</td><td>----</td><td>----</td><td>----</td><td>JH720449.1</td></tr>
<tr bgcolor="#EB8B31"><td>dd8730d9d33765ff135fcfadb8810280</td><td>----</td><td>----</td><td>----</td><td>GL383575.2</td></tr>
<tr bgcolor="#EB8B31"><td>df3e809f9a87f792218db18db51f6ad4</td><td>----</td><td>----</td><td>----</td><td>GL383522.1</td></tr>
<tr bgcolor="#EB8B31"><td>e0da36f2d1d2c6092f13d5bee52537e0</td><td>----</td><td>----</td><td>----</td><td>GL383516.1</td></tr>
<tr bgcolor="#EB8B31"><td>e0e934bd79ff323b31f4c9b80fb80a5c</td><td>----</td><td>----</td><td>----</td><td>JH159134.2</td></tr>
<tr bgcolor="#EB8B31"><td>e11adfbb638e60f61d7e8ef6647f30f2</td><td>----</td><td>----</td><td>----</td><td>JH720448.1</td></tr>
<tr bgcolor="#EB8B31"><td>e2cd68e2099fbd7cee557d6a7910768f</td><td>----</td><td>----</td><td>----</td><td>KB021646.2</td></tr>
<tr bgcolor="#EB8B31"><td>e363729ea23dad7c6802e7b439b4f668</td><td>----</td><td>----</td><td>----</td><td>KE332497.1</td></tr>
<tr bgcolor="#EB8B31"><td>e5b96eb9510763261839281c198607dd</td><td>----</td><td>----</td><td>----</td><td>GL582979.2</td></tr>
<tr bgcolor="#EB8B31"><td>e5cd94b0e0668debf81b82f405597b28</td><td>----</td><td>----</td><td>----</td><td>GL582970.1</td></tr>
<tr bgcolor="#EB8B31"><td>e6c232469067e8cadfa852a2ea5513b7</td><td>----</td><td>----</td><td>----</td><td>GL949749.1</td></tr>
<tr bgcolor="#EB8B31"><td>e8c870267b2a5261edb9d51d0efd6469</td><td>----</td><td>----</td><td>----</td><td>GL582975.1</td></tr>
<tr bgcolor="#EB8B31"><td>ebf72aeb4d53f0fd56e2e72967751f8a</td><td>----</td><td>----</td><td>----</td><td>GL383562.1</td></tr>
<tr bgcolor="#EB8B31"><td>ed6bcd4459b3bc6b366ce00262952f57</td><td>----</td><td>----</td><td>----</td><td>GL949743.1</td></tr>
<tr bgcolor="#EB8B31"><td>ed6fb45e0a25c31903cbb0f78d9d487e</td><td>----</td><td>----</td><td>----</td><td>GL383546.1</td></tr>
<tr bgcolor="#EB8B31"><td>edf086bce359065367b105cae0abfeee</td><td>----</td><td>----</td><td>----</td><td>GL383579.1</td></tr>
<tr bgcolor="#EB8B31"><td>edf41bfaf2584364bb4c5a645d73d53c</td><td>----</td><td>----</td><td>----</td><td>GL383553.2</td></tr>
<tr bgcolor="#EB8B31"><td>f2bfb99f84f2dd2ea538fe69ee786a0d</td><td>----</td><td>----</td><td>----</td><td>GL383525.1</td></tr>
<tr bgcolor="#EB8B31"><td>f486a5a44493d2e6bf72bf95ae898e3c</td><td>----</td><td>----</td><td>----</td><td>JH806574.2</td></tr>
<tr bgcolor="#EB8B31"><td>f7ee47af8d462cd9aeb6d40de99acb36</td><td>----</td><td>----</td><td>----</td><td>JH806585.1</td></tr>
<tr bgcolor="#EB8B31"><td>fa5fa49d281fc855dd1076c4f51bd8dc</td><td>----</td><td>----</td><td>----</td><td>KE332498.1</td></tr>
<tr bgcolor="#EB8B31"><td>faa48b73103366d1da02065870a58bda</td><td>----</td><td>----</td><td>----</td><td>GL383543.1</td></tr>
<tr bgcolor="#EB8B31"><td>faae4c952e9c38254538e1853b786276</td><td>----</td><td>----</td><td>----</td><td>GL582976.1</td></tr>
<tr bgcolor="#EB8B31"><td>fc93038463f9660e139435537ef53a5c</td><td>----</td><td>----</td><td>----</td><td>KB663604.1</td></tr>
<tr bgcolor="#EB8B31"><td>fdeb8db11e8544a638179a592c051331</td><td>----</td><td>----</td><td>----</td><td>JH636061.1</td></tr>
<tr bgcolor="#EB8B31"><td>fef0bc815f4826ea408515d8ec74ca80</td><td>----</td><td>----</td><td>----</td><td>JH636057.1</td></tr>
<tr bgcolor="#EB8B31"><td>ff7c4316cb69a8d571bd7ef85c1a10e4</td><td>----</td><td>----</td><td>----</td><td>JH720455.1</td></tr>
</table>

This table indicates that while most contigs contain the same data, there are several with **sequence differences between the references**.  Among those are **Chromosome 3**, **Chromosome Y**, and the **Mitochondrial Contig**.  

Anecdotally the changes are for bases for which there was low confidence, with those low confidence bases masked out to be the `IUPAC` symbol for any base.  However, _there does not seem to be a detailed comparison readily available (i.e. there's no proof that this is true)_.

Therefore, when doing comparisons across the four reference versions for each of these contigs, some care should be taken.

<a name="additional"></a>
## Additional Information
Some further details can be found on [this DNAnexus wiki page](https://wiki.dnanexus.com/scientific-notes/human-genome) as of 2019/01/30.

[Back to Top](#0) 

<hr>
