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
				gene_id = re.search(r">(\S+)(\s+)?", line) # procura sequencias que começam com '>', podendo ter outros caracteres sem ser espaços
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

all_frames = reading_frames(genes_dict) # Com isso, usamos a função e armazenamos sua saída no dicionário all_frames


# Criando a função para traduzir cada reading frame de cada gene
def codons_translation (all_frames): # neste bloco de for, os codons de cada frame, para cada gene são percorridos, e ocorre a formação da sequencia de proteina
	for gene_id in all_frames:
		protein_dict[gene_id] = {}
		for frames in all_frames[gene_id]:
			protein_dict[gene_id][frames] = []
			protein = ""
			for codon in all_frames[gene_id][frames]:
				if codon in translation_table:
					protein += translation_table[codon] #armazenando os aminoácidos traduzidos na string protein
				else:
					print (" O códon não se encontra no dicionário de tradução")
					exit()
			protein_dict[gene_id][frames] = protein

	for gene_id in protein_dict:
		longest_ORF[gene_id] = ""
		for frames in protein_dict[gene_id]:
			ORFs = re.finditer(r"(M[A-Z]+?\*)" , protein_dict[gene_id]) #.finditer() pesquisa todas as sequencias que começam com M e finalizam com * (códons de início e parada). Retorna um objeto com todas as ocorrências, que é iterável
			for i in ORFs:
				if len(i.group(1)) > len(longest_ORF[gene_id]):
					longest_ORF[gene_id] = i.group(1)
	return longest_ORF

longest_ORF = codons_translation(all_frames) # Com isso, usamos a função para traduzir os códons e armazenamos a saída no dicionário de proteínas

print (longest_ORF)

# Procurando pela maior sequencia de proteina em cada frame
# Considerando que o início é pela metionina ( representado por 'M' ) e que os códons de parada são representados por '*'

#for genes in protein_dict:
#	longest_protein[genes] = ""
#	proteinas = re.finditer(r"(M[A-Z]+?\*)" , protein_dict[genes]) # o .finditer() pesquisará todas as ocorrencias de sequencias começando com M, tendo uma ou mais letras maiúsculas e finalizando com *, na nossa string. Com isso, retornará um objeto iterável, que corresponde a todas as ocorrências encontradas, sendo que estas ocorrencias poderao ser percorridas por um loop for
#	for i in proteinas:
#		if len(i.group(1)) > len(longest_protein[genes]):
#			longest_protein[genes] = i.group(1)





# Escrevendo o resultado no arquivo de saída

#with open("teste_saida_08.txt" , "w") as output_file:  #uso de dois blocos de for
#	for gene_id , frames, in all_frames.items():  	#método items() retorna todos os items do dicionário como tuplas, neste caso transforma-os em um par (gene_id , frames), onde gene_id é a chave e 'frames' é o valor associado a essa chave (um dicionário que contem os frames 1,2,3 e seus respectivos códons). O for será usado para percorre pela nome dos genes e obter acesso ao dicionário frames
#		for frame, codons in frames.items():  #da mesma maneira que o anterior, aqui o método items() fará um par (frame, codons), onde frame será as chaves (1,2 ou 3) e 'codons' receberá os valores associados a cada chave, durante cada iteração. Ou seja, o for percorrerá cada frame, e atribuirá os valores aos códons
#			headline = f"{gene_id}-frame-{frame}-codons\n"
#			output_file.write (headline)
#			output_file.write (" ".join(codons) + "\n\n")


