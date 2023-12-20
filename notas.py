#!/usr/bin/env python3

# Script para calcular média de notas

# Inicializando variáveis
total = 0
contador_notas = 0

while contador_notas <= 9:
	contador_notas += 1 # Adicionando incremento para a variável
	try:
		nota = float(input(f"Informe a nota da Atividade {contador_notas}: ")) # Recebendo input do usuário
	except:
		print(f"Deve ser informado um número entre 0 e 10 para a nota da Atividade {contador_notas}!") # Forncecendo mensagem de erro
		exit()
	if nota <= 10 and nota >=0: # Condições para os valores limites da nota
		total += nota
	else:
		print(f"A nota da Atividade {contador_notas} deve ser um valor entre 0 e 10!") # Fornecendo mensagem de erro
		exit()

media_disciplina = total / 10 # Calculando media da disciplina

print(f"A média final da disciplina é: {media_disciplina}") # Fornecendo o valor de média como output ao usuário
