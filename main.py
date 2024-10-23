import pysam
# from Bio import SeqIO
# from Bio.Seq import Seq
from collections import defaultdict
import argparse
import os
from scripts.terminalPrinting import *
import pprint as pp
from csv import writer


def parseBam(bamfile, mappedReads, softClipReads, top_n = -1):
  bam = pysam.AlignmentFile(bamfile, "rb", threads = 20)
  readIndex = 0

  for read in bam:
    if readIndex % 10000 == 0:
      print("Parsing {}th read".format(str(readIndex)), end = "\r")

    if top_n != -1 and readIndex > top_n:
      print()
      return

    readIndex += 1

    # ignore if optical/PCR duplicate OR without a mate
    if (read.flag & 1024) or (not read.flag & 1):
      readIndex += 1
      continue

    # ignore if unmapped or mate is unmapped
    if (read.flag & 4) or (read.flag & 8):
      readIndex += 1
      continue
    
    refnameIsProviral = read.reference_name == "hivJLat10-6"
  
    cigarString = read.cigartuples
    # 4 is soft clip
    hasSoftClipAtEnd = cigarString != None and (cigarString[-1][0] == 4 or cigarString[0][0] == 4)

    # ignore if no valid barcode
    if not read.has_tag("CB"):
      continue

    if refnameIsProviral and not hasSoftClipAtEnd:
      mappedReads[read.query_name].append(read)
    elif refnameIsProviral and hasSoftClipAtEnd:
      softClipReads[read.query_name].append(read)

  print()
  return bam

def main(args):
  mappedReads = defaultdict(list)
  softClipReads = defaultdict(list)

  parseBam(args.bamfile, mappedReads=mappedReads, softClipReads=softClipReads, top_n = -1)

  # check where the reads are aligning to
  compiledData = []
  for pair in mappedReads:
    pairReal = mappedReads[pair]
    if len(pairReal) != 2:
      continue

    # needs valid cb to pass
    if pairReal[0].get_tag("CB")[-2:] != "-1":
      continue


    compiledData.append({
      "qname": pairReal[0].query_name,
      "cbc": pairReal[0].get_tag("CB"),
      "forward_refStart": pairReal[0].reference_start,
      "rev_refStart": pairReal[1].reference_start,
      "forward_length": pairReal[0].query_length,
      "rev_length": pairReal[1].query_length,
      "forward_cigar": pairReal[0].cigarstring,
      "rev_cigar": pairReal[1].cigarstring,
      "forward_alts": pairReal[0].get_tag("XA") if pairReal[0].has_tag("XA") else None,
      "rev_alts": pairReal[1].get_tag("XA") if pairReal[1].has_tag("XA") else None,
      "needToReorder": False
    })

  # note that these are 0-index
  # TODO use only the start position of the nested primer because a tagmented sample could
  # result in a shorter (or longer) segment on the 5' end.
  # could probably set up a buffer area to incorporate/filter
  HIV_NESTED_SETS = {
    "Set1": 747,
    "Set2": 2554, # double check this
    "Set3": 7794,
    "Set4": 10063
  }

  # BUFFER = 1
  
  alignedData = defaultdict(list)

  # note the 0-index correction from 1-index default in SAM/BAM files
  getPosFromAlt = lambda x: abs(int(x.split(",")[1])) - 1
  for d in compiledData:
    # use only reads with M (no insertion/deletion)
    # TODO will need to upgrade this for non cell line studies because
    # of mutation
    undesiredCigar = r"IDS"
    if any(ele in undesiredCigar for ele in d["forward_cigar"]):
      continue
    elif any(ele in undesiredCigar for ele in d["rev_cigar"]):
      continue

    primary = [d["forward_refStart"], d["rev_refStart"]]

    if d["forward_alts"] is not None and d["rev_alts"] is not None:
      secondary = [getPosFromAlt(d["forward_alts"]), getPosFromAlt(d["rev_alts"])]
    elif d["forward_alts"] is not None:
      secondary = [getPosFromAlt(d["forward_alts"]), d["rev_refStart"]]
    elif d["rev_alts"] is not None:
      secondary = [d["forward_refStart"], getPosFromAlt(d["rev_alts"])]
    else:
      secondary = None

    # reorder as needed because max alignment should be < 300
    maxAlignmentLen = 300
    if abs(primary[1] - primary[0]) > maxAlignmentLen and secondary is not None and secondary[1] - primary[0] < maxAlignmentLen :
      old = primary[1]
      primary[1] = secondary[1]
      secondary[1] = old
      d["needToReorder"] = True

    addToggle = None
    for hsets in HIV_NESTED_SETS:
      if primary[1] == HIV_NESTED_SETS[hsets] or primary[1] == HIV_NESTED_SETS[hsets] - 1:
        alignedData[d["cbc"]].append(hsets)
        addToggle = "P"

      elif secondary is not None and (secondary[1] == HIV_NESTED_SETS[hsets] or secondary[1] == HIV_NESTED_SETS[hsets] - 1):
        alignedData[d["cbc"]].append(hsets)
        addToggle = "S"

      if addToggle is not None:
        break
    

    if args.debug:
      if addToggle is None:
        printRed((d["qname"], primary, secondary, addToggle, d["needToReorder"], hsets))
      else:
        printGreen((d["qname"], primary, secondary, addToggle, d["needToReorder"], hsets))

  # export
  outputFNs = {
    "cbc_sets": "cbc_sets.tsv"
  }

  for k in outputFNs:
    outputFNs[k] = args.outputDir + "/" + outputFNs[k]

  with open(outputFNs["cbc_sets"], "w") as tsvfile2:
    writ2 = writer(tsvfile2, delimiter = "\t")

    writ2.writerow(["cbc", "sets"])
    for data in alignedData:
      sets = ",".join(set(alignedData[data]))
      writ2.writerow([data, sets])


if __name__ == '__main__':
  # set up command line arguments
  parser = argparse.ArgumentParser(
    description = "Identify HIV-associated reads amidst cellular scATAC reads")

  parser.add_argument("--bamfile",
    required = True,
    help = "Aligned and namesorted BAM file for GOTCHAseq")
  parser.add_argument("--outputDir",
    required = True,
    help = "Output file location")
  parser.add_argument("--debug",
    action = 'store_const',
    default = False,
    const = True,
    help = "Debugging mode; high verbosity")

  args = parser.parse_args()

  if not os.path.exists(args.outputDir):
    os.makedirs(args.outputDir)

  if not os.path.exists(args.bamfile):
    raise Exception("BAM file not found")

  main(args)

