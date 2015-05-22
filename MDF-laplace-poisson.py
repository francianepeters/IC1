#!/usr/bin/env python
# *-* coding: utf-8 *-*
#Teste do meld
import numpy 
from scipy import linalg
import matplotlib.pyplot as plt

# Franciane: coloquei esse comentario so para editar o codigo e enviar uma segunda versao para o repositorio


#Imprime uma lista em formato de matriz
def imprimeMatriz(matriz):
    print'\n'
    for i in range (len(matriz)):
        for j in range (len(matriz[0])):
            print matriz[i][j],
        print

#Imprime uma lista em formato de vetor
def imprimeVetor(vetor):
    print'\n'
    for i in range (len(vetor)):
        print vetor[i]

#Avalia se o problema pode ser resolvido
#Tc1 e Tc2 nao podem ser naturais simultaneamente
def verificacc(tc1,tc2):
    
    if tc1 != 1 and tc1 != 0:
        print 'Erro cc 1'
        
    if tc2 != 1 and tc2 != 0:
        print 'Erro cc 2'
        
    if tc1 == 1 and tc2 == 1:
        print 'Ambas as condicoes de contorno sao naturais => impossivel definir solucao'
        
    
#Avalia a entrada do numero de pontos
#Deve ser um valor inteiro >= 2
def verificanpts(npontos):
    if isinstance(npontos,int)== False:
        print 'Erro no numero de pontos'
    if npontos < 2:
        print 'Erro no numero de pontos.'

#Montar a matriz de coeficientes da equação da laplace
#Via MDF
def taylorcoef1d(npontos,tc1,vc1,tc2,vc2):  
    h=1.0/(npontos-1) #h-dimensao dos trechos particionados

    matriz=[0]*npontos #linhas da matriz->cria uma lista com (npontos) elementos

    for i in range (npontos):
        matriz[i]=[0]*npontos #colunas da matriz->substitui cada elemento da lista por uma lista com (npontos) elementos

    if tc1==0: # 1a linha: condição de contorno essencial do ponto inicial       
        matriz[0][0]=1.0 
    elif tc1==1:# 1a linha: condição de contorno natural do ponto inicial
        matriz[0][0]=-1.0/h
        matriz[0][1]= 1.0/h


#preenche a linha 1 e a linha [(npontos)-1] da matriz de coeficientes  
    for i in range(1,npontos-1):
        matriz[i][i+1]=1/h**2
        matriz[i][i-1]=1/h**2
        matriz[i][i]=-2/h**2

    if tc2==0: # última linha: condição de contorno essencial do ponto final       
        matriz[npontos-1][npontos-1]=1.0
    elif tc2==1:# ultima linha: condição de contorno natural do ponto final
        matriz[npontos-1][npontos-2]=-1.0/h
        matriz[npontos-1][npontos-1]= 1.0/h

    return matriz

#Monta o vetor C da equacao matricial (taylorcoef1d)X = C
#Sendo X o vetor com os valores de x onde se deseja definir u(x)
def taylorvetor(npontos,tc1,vc1,tc2,vc2):
    vetor=[0]*npontos

    vetor[0]=vc1
    vetor[npontos-1]=vc2

    return vetor

#Montar o vetor resposta U(x) para os pontos x desejados
#Laplace solucionado via analitica
def exata1d(npontos,tc1,vc1,tc2,vc2):  
    h=1.0/(npontos-1)
    exata=[0]*npontos
       
    if tc1==0 and tc2==0:
         for i in range (npontos):
            x=i*h
            exata[i]=(vc2-vc1)*x+vc1

    if tc1==0 and tc2==1:
         for i in range (npontos):
            x=i*h
            exata[i]=vc2*x+vc1

    if tc1==1 and tc2==0:
         for i in range (npontos):
            x=i*h
            exata[i]= vc1*x+(vc2-vc1)

    return exata





    
npontos=input('número de pontos para discretizacao:')

verificanpts(npontos)

tc1=input('tipo de condicao de contorno da esquerda:')
vc1=input('valor de condicao de contorno da esquerda:')

tc2=input('tipo de condicao de contorno da direita:')
vc2=input('valor de condicao de contorno da direita:')

verificacc(tc1,tc2)

matrizcoeftaylor=taylorcoef1d(npontos,tc1,vc1,tc2,vc2)

imprimeMatriz(matrizcoeftaylor)
print 'Matriz de Coeficientes Taylor'

vetortaylor=taylorvetor(npontos,tc1,vc1,tc2,vc2)

imprimeVetor(vetortaylor)
print 'Vetor Taylor'

solnumerica=linalg.solve(matrizcoeftaylor,vetortaylor)

imprimeVetor(solnumerica)
print 'Soluçao Numerica'

abscissas=numpy.linspace(0,1,npontos)

imprimeVetor(abscissas)
print 'Pontos discretizados'

solexata=exata1d(npontos,tc1,vc1,tc2,vc2)

imprimeVetor(solexata)
print 'Soluçao Exata'


plt.plot(abscissas,solnumerica,'bo')   # estou a fazer o gráfico dos dois vectores.

    
plt.plot(abscissas,solexata,'rx')   # estou a fazer o gráfico dos dois vectores.


plt.show()    # Esta linha de código serve para mostrar o gráfico total.
