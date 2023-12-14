#!/usr/bin/env python3

# Script para calcular média de notas

# Inicializando variáveis
total = 0
contador_notas = 0

while contador_notas <= 9:
	contador_notas += 1
	try:
		nota = float(input(f"Informe a nota da Atividade {contador_notas}: "))
	except:
		print(f"Deve ser informado um número entre 0 e 10 para a nota da Atividade {contador_notas}!")
		exit()
	if nota <= 10 and nota >=0:
		total += nota
	else:
		print(f"A nota da Atividade {contador_notas} deve ser um valor entre 0 e 10!")
		exit()

media_disciplina = total / 10

print(f"A média final da disciplina é: {media_disciplina}")
