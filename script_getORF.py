#!/usr/bin/env python3

# Importando bibliotecas
import sys
import re

# Inserindo tabela de tradução de códons
translation_table = {
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
    'AAT':'N', 'AAC':'N',
    'GAT':'D', 'GAC':'D',
    'TGT':'C', 'TGC':'C',
    'CAA':'Q', 'CAG':'Q',
    'GAA':'E', 'GAG':'E',
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
    'CAT':'H', 'CAC':'H',
    'ATT':'I', 'ATC':'I', 'ATA':'I',
    'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    'AAA':'K', 'AAG':'K',
    'ATG':'M',
    'TTT':'F', 'TTC':'F',
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'TGG':'W',
    'TAT':'Y', 'TAC':'Y',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    'TAA':'*', 'TGA':'*', 'TAG':'*'
}

# Inicializando variáveis
genes_dict = {}
all_frames = {}
genes_rev_dict = {}
longest_ORF = {}
protein_dict = {}

# Recebendo o caminho do arquivo multi-FASTA que será lido
try:
	file_path = sys.argv[1]
except IndexError:
	print("Forneça o nome de um arquivo após ./script_getORF.py") # mensagem de erro para quando o usuário não fornece um arquivo de input

# Realizando a leitura do arquivo multi-FASTA
try:
	with open(file_path, "r") as file:
		for line in file:
			if line.startswith('>'):
				line = line.rstrip() # remove quebra de linha
				gene_id = re.search(r">(\S+)(\s+)?", line) # procura sequencias que começam com '>', seguido de outros caracteres sem espaços em branco
				gene_id = gene_id.group(1) # armazena o nome do gene
				genes_dict[gene_id] = "" # cria chaves do dicionário
				full_gene_seq = "" # cria a string vazia, que armazenara sequencia do gene. A cada loop, ela é reinicializada
			else:
				full_gene_seq += line.rstrip().upper() # remove quebra de linha, deixa nucleotideos em letra maiúscula e junta as linhas em uma string só
				genes_dict[gene_id] = full_gene_seq # armazena a sequencia completa de cada gene
except IOError:
	print("Não foi possível encontrar o arquivo:" , file_path) # mensagem de erro para quando o arquivo de input não é encontrado

# Criando a função para calcular as 6 fases de leitura
def reading_frames (all_genes):
	for genes in all_genes:
		all_frames[genes] = {'1' : [] , '2' : [] , '3' : [] , '4' : [] , '5' : [] , '6' : [] } # criando as chaves no dicionário, as quais receberão a lista de códons de cada fase de leitura
		all_frames[genes]['1'] = re.findall(r"(.{3})", all_genes[genes]) # cálculo da fase de leitura 1
		all_frames[genes]['2'] = re.findall(r"(.{3})", all_genes[genes][1:]) # cálculo da fase de leitura 2
		all_frames[genes]['3'] = re.findall(r"(.{3})", all_genes[genes][2:]) # cálculo da fase de leitura 3
		genes_rev_dict[genes] = all_genes[genes][::-1] # armazena o reverso da sequencia de cada gene
		genes_rev_dict[genes] = genes_rev_dict[genes].replace('G','c').replace('C','g').replace('T','a').replace('A','t').upper() # armazena o reverse complement de cada gene
		all_frames[genes]['4'] = re.findall(r"(.{3})", genes_rev_dict[genes]) # cálculo da fase de leitura 4
		all_frames[genes]['5'] = re.findall(r"(.{3})", genes_rev_dict[genes][1:]) # cálculo da fase de leitura 5
		all_frames[genes]['6'] = re.findall(r"(.{3})", genes_rev_dict[genes][2:]) # cálculo da fase de leitura 6
	return all_frames

all_frames = reading_frames(genes_dict) # Com isso, usamos a função acima e armazenamos sua saída no dicionário all_frames

# Criando a função para fornecer o ORF mais comprido de cada gene, bem como sua posição de inicio e final, e qual frame está localizado
def long_ORF (all_frames):
	for gene_id in all_frames:
		longest_ORF[gene_id] = {"sequence" : "", "start" : 0 , "end" : 0 , "frame" : 0, "len" : 0} # cria um dicionário com as informações que precisarão ser armazenadas para cada gene
		for frames in all_frames[gene_id]:
			sequence = ''.join(all_frames[gene_id][frames]) # junta a lista de códons de cada frame em uma string só, para que possamos usar .finditer. Não pode separar os codons por espaços para que tenhamos as posições de inicio e final corretas
			ORFs = re.finditer(r"(ATG[A-Z]*?(TAA|TGA|TAG))", sequence) # Pesquisa todas as sequencias de DNA que começam com ATG e terminam com algum dos códons de parada. Armazena as correspondências em um objeto iterável
			for i in ORFs: # usado para percorrer cada correspondência encontrada
				orf_sequence = i.group() # armazena cada correspondência
				if len(orf_sequence) > len(longest_ORF[gene_id]["sequence"]) and len(orf_sequence) % 3 == 0: # caso o ORF atual seja maior que o anterior e também seja uma sequência múltipla de 3, ele é armazenado, junto com sua posição de inicio e final de determinado frame
					longest_ORF[gene_id]["sequence"] = orf_sequence
					longest_ORF[gene_id]["start"] = i.start() +1 # soma um para ter a posição na sequencia de DNA
					longest_ORF[gene_id]["end"] = i.end() +1
					longest_ORF[gene_id]["frame"] = int(frames)
					longest_ORF[gene_id]["len"] = len(orf_sequence)
	return longest_ORF

longest_ORF = long_ORF(all_frames) # Com isso, usamos a função acima e armazenamos sua saída no dicionário longest_ORF

# Criando a função para traduzir os ORFs mais compridos de cada gene
def codons_translation (longest_ORF):
	for gene_id in longest_ORF:
		protein_dict[gene_id] = {"sequence" : "", "start" : longest_ORF[gene_id]["start"] , "end" : longest_ORF[gene_id]["end"] , "frame" : longest_ORF[gene_id]["frame"]}
		codons = []
		protein = ""
		codons = re.findall(r"(.{3})", longest_ORF[gene_id]["sequence"]) # divide a sequencia do ORF em códons
		for codon in codons:
			if codon in translation_table:
				protein += translation_table[codon] # armazenando os aminoácidos traduzidos na variável protein (string)
			else:
				print (f"O códon {codon} não se encontra no dicionário de tradução") # mensagem de erro caso um códon não seja contemplado pela lista de tradução
				sys.exit()
			protein_dict[gene_id]["sequence"] = protein
	return protein_dict

protein_dict = codons_translation(longest_ORF) # Com isso, usamos a função acima e armazenamos sua saída no dicionário protein_dict

# Escrevendo os ORFs mais compridos no arquivo de saída
with open("ORF.fna" , "w") as output_file:
	for gene_id in longest_ORF:
		frame = longest_ORF[gene_id]["frame"] # Definindo as variáveis que serão escritas
		start = longest_ORF[gene_id]["start"]
		end = longest_ORF[gene_id]["end"]
		sequence = longest_ORF[gene_id]["sequence"]

		headline = f">{gene_id}_frame{frame}_{start}_{end}\n" # Escreve o cabeçalho no formato FASTA iniciado por >
		output_file.write (headline)

		for i in range(0, len(sequence), 60): # Escreve linhas de sequencia de DNA com no maximo 60 caracteres
		    output_file.write(sequence[i:i+60] + "\n")

# Escrevendo os peptideos traduzidos de cada ORF
with open("ORF.faa" , "w") as output_file:
        for gene_id in protein_dict:
                frame = protein_dict[gene_id]["frame"] # Definindo as variáveis que serão escritas
                start = protein_dict[gene_id]["start"]
                end = protein_dict[gene_id]["end"]
                sequence = protein_dict[gene_id]["sequence"]

                headline = f">{gene_id}_frame{frame}_{start}_{end}\n" # Escreve o cabeçalho no formato FASTA iniciado por >
                output_file.write (headline)

                for i in range(0, len(sequence), 60): # Escreve linhas de sequencia de DNA com no maximo 60 caracteres
                    output_file.write(sequence[i:i+60] + "\n")

