shell.prefix("source $HOME/.bashrc; ")
shell.prefix("source /etc/profile; ")

IDS, = glob_wildcards("genomes/{id}_genomic.fna.gz")

configfile: "config.json"

ruleorder: ncbi2table > gff2table > proteins

rule all:
        input: expand('puls/{sample}.puls.sum.tsv', sample=IDS), expand('puls/{sample}.puls.tsv', sample=IDS), expand('dbcan/{sample}.out.dm.ps.filtered', sample=IDS), expand('pfam/{sample}.pfam', sample=IDS)

rule proteins:
	input: 'genomes/{id}_genomic.fna.gz'
	output: 
		faa='proteins/{id}_protein.faa',
		gff='proteins/{id}_prodigal.gff'
	conda: "envs/PULpy.yaml"
	shell: 'zcat {input} | prodigal -a {output.faa} -q -f gff > {output.gff}'

rule gff2table:
	input: "proteins/{id}_protein.faa", "proteins/{id}_prodigal.gff"
	output: 'feature_table/{id}_ft.txt'
	run:

		# if output dir doesn't exist
		# create it
		try:
    			os.stat("feature_table")
		except:
    			os.mkdir("feature_table")

		# open the input file
		gff = open(input[1], mode="r")

		# open the output file
		out = open(output[0], mode="w")

		# iterate over file
		for row in gff:

			if row.startswith("#"):
				continue

			arr = row.split('\t')
			if (len(arr) < 9):
				continue

			if (arr[2] != "CDS"):
				continue

			lstr = arr[8]
			lste = lstr.split(";")
			pnum = lste[0].split("_")[1]

			out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (arr[0],arr[3],arr[4],arr[6],arr[0]+"_"+pnum,arr[0]+"_"+pnum))

		gff.close();
		out.close();


rule ncbi2table:
	input: "ncbi_feature_table/{id}_feature_table.txt", "proteins/{id}_protein.faa", "proteins/{id}_genomic.gff"
	output: 'feature_table/{id}_ft.txt'
        run:

		# if output dir doesn't exist
		# create it
		try:
			os.stat("feature_table")
		except:
			os.mkdir("feature_table")

		# open the input file
		gff = open(input[0], mode="r")

		# open the output file
		out = open(output[0], mode="w")

		# skip one lines
		row1 = gff.readline()

		# iterate over file
		for row in gff:

			arr = row.split('\t')
			if (len(arr) < 20):
				continue

			if (arr[0] != "CDS"):
				continue

			out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (arr[6],arr[7],arr[8],arr[9],arr[10],arr[16]))

		gff.close();
		out.close();


rule pfam_scan:
	input: 'proteins/{id}_protein.faa'
	output:'pfam/{id}.pfam'
	threads: 4
	conda: "envs/PULpy.yaml"
	params:
		pfam=config["pfam_dir"]
	shell: "pfam_scan.pl -outfile {output} -as -cpu {threads} -fasta {input} -dir {params.pfam}"

rule dbcan:
	input: 'proteins/{id}_protein.faa'
	output:
		dm='dbcan/{id}.out.dm',
		out='dbcan/{id}.out'
	conda: "envs/PULpy.yaml"
	params:
		hmm=config["dbcan_hmm"]
	shell: "hmmscan --domtblout {output.dm} {params.hmm} {input} > {output.out}"

rule dbcan_filter:
	input: 'dbcan/{id}.out.dm'
	output:
		ps='dbcan/{id}.out.dm.ps',
		filt='dbcan/{id}.out.dm.ps.filtered'
	params:
		script=config["dbcan_dir"]+"/"+"hmmscan-parser.sh"
	shell: "{params.script} {input} > {output.ps} && cat {output.ps} | awk '$10 >= 0.35' | awk '$5 <= 1e-18' > {output.filt}"

rule puls:
	input: 
		pfam='pfam/{id}.pfam',
		dbcan='dbcan/{id}.out.dm.ps.filtered',
		ft='feature_table/{id}_ft.txt'
	params:
		id="{id}"
	output:
		all='puls/{id}.puls.tsv',
		sum='puls/{id}.puls.sum.tsv'
	conda: "envs/PULpy.yaml"
	shell: 
		'''
		mkdir -p puls
		./scripts/predict_puls.R {input.ft} {input.pfam} {input.dbcan} {output.all} {output.sum} {params.id}
		'''

