#!/usr/bin/env python
# *-* coding: utf-8 *-*

import numpy 
from scipy import linalg
import matplotlib.pyplot as plt


#Imprimi uma lista em formato de matriz
def imprimiMatriz(matriz):
    print'\n'
    for i in range (len(matriz)):
        for j in range (len(matriz[0])):
            print matriz[i][j],
        print
    print'\n'

#Imprime uma lista em formato de vetor
def imprimiVetor(vetor):
    print'\n'
    for i in range (len(vetor)):
        print vetor[i]

        
#Montar a matriz da equação da laplace
def laplace1d(npontos,tc1,vc1,tc2,vc2):  
    h=1.0/(npontos-1)

    matriz=[0]*npontos
    for i in range (npontos):
        matriz[i]=[0]*npontos

    if tc1==0: # 1a linha: condição de contorno essencial do ponto inicial       
        matriz[0][0]=1.0 
    elif tc1==1:# 1a linha: condição de contorno natural do ponto inicial
        matriz[0][0]=-1.0/h
        matriz[0][1]= 1.0/h
    else:
        print'ERRO tipo cc inicial'

        


      
    for i in range(1,npontos-1):
        matriz[i][i+1]=1/h**2
        matriz[i][i-1]=1/h**2
        matriz[i][i]=-2/h**2


    if tc2==0: # última linha: condição de contorno essencial do ponto final       
        matriz[npontos-1][npontos-1]=1.0
    elif tc2==1:# ultima linha: condição de contorno natural do ponto final
        matriz[npontos-1][npontos-2]=-1.0/h
        matriz[npontos-1][npontos-1]= 1.0/h
    else:
        print'ERRO tipo CC final'
     
        
    # última linha: condição de contorno do ponto final
    matriz[npontos-1][npontos-1]=1.0	 

    return matriz

#Montar a matriz da equação da laplace
def exatalaplace1d(npontos,tc1,vc1,tc2,vc2):  
    h=1.0/(npontos-1)
    exata=[0]*npontos
       
    if tc1==0 and tc2==0:
         for i in range (npontos):
            x=i*h
            exata[i]=(vc2-vc1)*x+vc1

    return exata


#Montar um vetor independente do sistema
def montavetor(npontos,tc1,vc1,tc2,vc2):
    vetor=[0]*npontos

    vetor[0]=vc1
    vetor[npontos-1]=vc2

    return vetor


    
npontos=input('número de pontos para discretizacao:')

#if(npontos<=1)then
#('número de pontos incorreto')

tc1=input('tipo de condicao de contorno da esquerda:')
vc1=input('valor de condicao de contorno da esquerda:')


tc2=input('tipo de condicao de contorno da direita:')
vc2=input('valor de condicao de contorno da direita:')


matriz=laplace1d(npontos,tc1,vc1,tc2,vc2)

imprimiMatriz(matriz)

vetor=montavetor(npontos,tc1,vc1,tc2,vc2)

imprimiVetor(vetor)

solucao=linalg.solve(matriz,vetor)

imprimiVetor(solucao)

x=numpy.linspace(0,1,npontos)

imprimiVetor(x)

exata=exatalaplace1d(npontos,tc1,vc1,tc2,vc2)

imprimiVetor(exata)


plt.plot(x,solucao,'bo')   # estou a fazer o gráfico dos dois vectores.

    
plt.plot(x,exata,'rx')   # estou a fazer o gráfico dos dois vectores.


plt.show()    # Esta linha de código serve para mostrar o gráfico total.
