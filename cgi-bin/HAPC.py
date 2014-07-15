#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin/python
import cgi, sys, math
from decimal import Decimal

form = cgi.FieldStorage()
hlas = form.getvalue("hlasname")
seqs = form.getvalue("sequencesname")
runHAPC = form.getvalue("runHAPC")

def printHtmlHeaders():
	print "Content-Type: text/html"
	print
	print """<!DOCTYPE html><html><head>
	<script src="//ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"></script>
	<script src="../cgi-bin/script.js"></script>
	<link rel="stylesheet" href="../css/style.css"></head><body>"""

def printFileHeaders(filename):
	print "Content-Disposition: attachment; filename=\""+filename+"\""
	print "Content-Type:application/octet-stream; name=\""+filename+"\""
	print

codon_dict = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
			  'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
			  'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
			  'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
			  'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
			  'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
			  'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
			  'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
			  'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
			  'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
			  'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
			  'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
			  'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
			  'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
			  'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
			  'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
			  '---':'-', 'XXX':'-', '???':'?'}

mixture_dict = {'W':'AT', 'R':'AG', 'K':'GT', 'Y':'CT',
				'S':'CG', 'M':'AC', 'V':'AGC', 'H':'ATC',
				'D':'ATG', 'B':'TGC', 'N':'ATGC', '-':'-'}

def resolveCodon(codon):
	nonmix = []
	if (codon in codon_dict):
		return [codon]
	elif (codon.count('-') + codon.count('X') == 3):
		return ['---']
	elif (1 <= codon.count('-') <= 2) or (1 <= codon.count('X') <= 2):
		return ['???']
	for base in codon:
		# Check for mixtures
		if (base in mixture_dict):
			if (not nonmix):
				nonmix = [x for x in mixture_dict[base]]
			else:
				nonmix = [x+y for x in nonmix for y in mixture_dict[base]]
		else:
			if (not nonmix):
				nonmix.append(base)
			else:
				nonmix = [x+base for x in nonmix]
	return nonmix

# Flag can be 0, 1, or 2 depending on the desired output
# Flag = 1 will output all mixtures as "X"
# Flag = 2 will output all synonymous mixtures as they are and all non-synonymous mixtures as "X"
# Flag = 3 will output all mixtures in the format [A/B] if a mixture encodes for amino acid A or B
def translateDNA(sequence, resolvecharacter="X", flag=2):
	sequence = sequence.translate(None, ' \n\r\n').upper()
	aaseq = []
	i = 0
	while i < len(sequence):
		codon = resolveCodon(sequence[i:i+3])
		# If the codon has no mixture bases just add it to the amino acid chain
		if len(codon) <= 1:
			aaseq.append(codon_dict[codon[0]])
		# Codon contains mixture base
		else:
			# If flag is set to 1
			if (flag == 1):
				aaseq.append(resolvecharacter)
			# If flag is set to 2
			elif (flag == 2):
				unique = set([codon_dict[potential] for potential in codon])
				# If there is more than resolved one amino acid
				if (len(unique) > 1):
					aaseq.append(resolvecharacter)
				else:
					aaseq.append(unique.pop())
			# If flag is set to 3
			else:
				unique = set([codon_dict[potential] for potential in codon])
				# If there is more than resolved one amino acid
				if (len(unique) > 1):
					aaseq.append('['+('/').join(unique)+']')
				else:
					aaseq.append(unique.pop())
		i += 3
	return aaseq

def parse(input):
	val = input.translate(None, '\r')
	val = [x.split('\t') for x in val.split('\n')]
	return val

def parseHLA(hla, res=4):
	rval = hla.strip()
	rval = hla.translate(None, ":*")
	try:
		int(rval[-1])
	except (ValueError, IndexError) as e:
		rval = rval[:-1]
	return rval[:res+1]

def groupHLA(hlas):
	rdic = {}
	for pair in hlas:
		hla = parseHLA(pair[0])
		loc = pair[1][:-1]
		aa = pair[1][-1]
		if hla not in rdic:
			rdic[hla] = [[loc, aa]]
		else:
			rdic[hla].append([loc, aa])
	return rdic

def getSeqs(seqs):
	d = {}
	for patient in seqs:
		pid = patient[0]
		nuseq = patient[-1]
		if (len(nuseq) % 3 != 0):
			aaseq = "Not divisible by 3"
		else:
			aaseq = translateDNA(nuseq)
		d[pid] = {'A': [], 'B': [], 'C': [], 'seq': aaseq}
		for hla in patient[1:-1]:
			hla = parseHLA(hla)
			if (hla == ""):
				continue
			if (hla[0].upper() == 'A'):
				d[pid]['A'].append(hla)
			elif (hla[0].upper() == 'B'):
				d[pid]['B'].append(hla)
			else:
				d[pid]['C'].append(hla)
	return d
	

def checkSeqs(interest, d, N):
	hla_and_mut, only_hla, only_mut, neither, line = 0,0,0,0,0
	hla = parseHLA(interest[0])
	pos = int(interest[1][:-1])
	mut = interest[1][-1]
	for patient in d.keys():
		# If there is no HLA data for this HLA type (A,B, or C)
		if (len(d[patient][hla[0].upper()]) == 0):
			return ["No data for hla {}".format(hla[0].upper())]*5
		else:
			pathlas = d[patient][hla[0].upper()]
		# Try to look at the position in the translated sequence, making sure not to go out of range
		try:
			comparebase = d[patient]['seq'][pos-1]
		except IndexError:
			return "Error: position out of range in sequence in patient {}".format(patient)
		# There are three situations where a patient will not be included in the analysis
		#
		# 1. The amino acid at the position of interest is unknown
		# 2. The patient is missing one or both HLA type[s] of interest (A,B, or C)
		#    and the HLA of interest is not present within the patient
		# 3. The HLA type of interest is a subtype of the HLA type of the patient
		if (comparebase not in codon_dict.values()):
			N -= 1
			continue
		minlen = min(len(hla), [len(x) for x in pathlas])
		if (hla in pathlas):
			if (mut == comparebase):
				hla_and_mut += 1
			else:
				only_hla += 1
		elif (hla[:minlen] in [x[:minlen] for x in pathlas]):
			if (len(hla) > len(max(pathlas))):
				N -= 1
				continue
			else:
				if (mut == comparebase):
					hla_and_mut += 1
				else:
					only_hla += 1
		else:
			if (len(pathlas) < 2):
				N -= 1
				continue
			if (mut == comparebase):
				only_mut += 1
			else:
				neither += 1
	return [hla_and_mut, only_hla, only_mut, neither, N]

if (runHAPC is not None):
	printHtmlHeaders()
	hlas = parse(hlas)
	seqs = parse(seqs)
	#print hlas,seqs
	d = getSeqs(seqs)
	N = len(d.keys())
	print """<div class="container"><table id="output_table">
	<tr>
		<td>Position</td>
		<td>Amino Acid</td>
		<td>HLA</td>
		<td>HLA+ AA+</td>
		<td>HLA+ AA-</td>
		<td>HLA- AA+</td>
		<td>HLA- AA-</td>
		<td>total</td>
		<td>lnOR</td>
		<td>p</td>
	</tr>
	"""
	for interest in hlas:
		hla = parseHLA(interest[0])
		pos = interest[1][:-1]
		mut = interest[1][-1]
		count = checkSeqs(interest, d, N)
		try:
			abfact = math.factorial(count[0]+count[1])
			cdfact = math.factorial(count[2]+count[3])
			acfact = math.factorial(count[0]+count[2])
			bdfact = math.factorial(count[1]+count[3])
			afact = math.factorial(count[0])
			bfact = math.factorial(count[1])
			cfact = math.factorial(count[2])
			dfact = math.factorial(count[3])
			nfact = math.factorial(count[4])
			logp = math.log(abfact) + math.log(cdfact) + math.log(acfact) + math.log(bdfact) - math.log((afact * bfact * cfact * dfact * nfact))
			pval = round(math.exp(logp),8)
		except ZeroDivisionError:
			pval = "Indeterminate"
		try:
			lnOR = round(math.log(float(count[0]*count[3])/(count[2]*count[3])),8)
		except (ZeroDivisionError, ValueError) as e:
			lnOR = "Indeterminate"
		print """
		<tr>
			<td>{}</td>
			<td>{}</td>
			<td>{}</td>
			<td>{}</td>
			<td>{}</td>
			<td>{}</td>
			<td>{}</td>
			<td>{}</td>
			<td>{}</td>
			<td>{:.8f}</td>
		</tr>
		""".format(pos, mut, hla, count[0], count[1], count[2], count[3], count[4], lnOR, pval)
	print "</table></div></body></html>"
