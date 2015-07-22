#!/usr/bin/env python
# *-* coding: utf-8 *-*
#https://www.wolframalpha.com/input/?i=u%3D%28%28%28x^4%29%2F12%29-%28%28x^3%29%2F6%29%29%2B%28%281%2F12%29%29x
#http://interactivepython.org/NWUbZ/LMdYZ/LpOMZ/courselib/static/pythonds/Graphs/graphintro.html

import numpy 
from scipy import linalg
import matplotlib.pyplot as plt



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
        print 'Erro no numero de pontos (nao inteiro)'
    if npontos < 2:
        print 'Erro no numero de pontos (<2)'

#Montar a matriz de coeficientes da equação de laplace
#Via MDF
def MDFcoef1d(npontos,tc1,vc1,tc2,vc2):  
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

#Monta o vetor C da equacao matricial (MDFcoef1d)X = C
#Sendo X o vetor com os valores de x onde se deseja definir u(x):

#Problema de Laplace
def vetorlaplace(npontos,tc1,vc1,tc2,vc2):
    vetor=[0]*npontos

    vetor[0]=vc1
    vetor[npontos-1]=vc2

    return vetor

#Problema de Poisson para f(x) = x² - x
def vetorpoisson(npontos,tc1,vc1,tc2,vc2):
    vetor=[0]*npontos
    h=1.0/(npontos-1)
    vetor[0]=vc1
    
    for i in range (1,npontos-2):
        x=i*h
        funcao=(x**2)-x
        vetor[i]=funcao
        
    vetor[npontos-1]=vc2

    return vetor

#Montar o vetor resposta U(x) para os pontos x desejados
#Laplace solucionado via analitica
def exatalaplace1d(npontos,tc1,vc1,tc2,vc2):  
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

#Montar o vetor resposta U(x) para os pontos x desejados
#Poisson solucionado via analitica para f(x) = x² - x
def exatapoisson1d(npontos,tc1,vc1,tc2,vc2):  
    h=1.0/(npontos-1)
    exata=[0]*npontos
       
    if tc1==0 and tc2==0:
        for i in range (npontos):
            x=i*h
            u=(((x**4)/12.0)-((x**3)/6.0))+((1.0/12.0)+vc2-vc1)*x + vc1
            exata[i]=u
            
    if tc1==0 and tc2==1:
        for i in range (npontos):
            x=i*h
            u=(((x**4)/12.0)-((x**3)/6.0))+(vc2 + (1/6.0))*x + vc1
            exata[i]=u

    if tc1==1 and tc2==0:
        for i in range (npontos):
            x=i*h
            u=(((x**4)/12.0)-((x**3)/6.0))+(vc1)*x + ((1/12.0)+vc2-vc1)
            exata[i]=u

    return exata

def erroabsoluto(vetorMDF,solexata,npontos):
    erro=[0]*npontos
    for i in range (npontos):
        erro[i]=vetorMDF[i]-solexata[i]
    return erro

##############################################################################
#Estrutura principal do programa:
##############################################################################

k=input('Digite "1" para resolver Laplace; "2" para Poisson (f(x)=x2-x):')
    
npontos=input('número de pontos para discretizacao:')

verificanpts(npontos)

tc1=input('tipo de condicao de contorno da esquerda:')
vc1=input('valor de condicao de contorno da esquerda:')

tc2=input('tipo de condicao de contorno da direita:')
vc2=input('valor de condicao de contorno da direita:')

verificacc(tc1,tc2)

matrizcoef=MDFcoef1d(npontos,tc1,vc1,tc2,vc2)

imprimeMatriz(matrizcoef)
print 'Matriz de Coeficientes MDF'

if k==1:
    vetorMDF=vetorlaplace(npontos,tc1,vc1,tc2,vc2)

elif k==2:
    vetorMDF=vetorpoisson(npontos,tc1,vc1,tc2,vc2)

imprimeVetor(vetorMDF)
print 'Vetor MDF'

solnumerica=linalg.solve(matrizcoef,vetorMDF)

imprimeVetor(solnumerica)
print 'Soluçao Numerica'

abscissas=numpy.linspace(0,1,npontos)

imprimeVetor(abscissas)
print 'Pontos discretizados'

if k==1:
    solexata=exatalaplace1d(npontos,tc1,vc1,tc2,vc2)

elif k==2:
    solexata=exatapoisson1d(npontos,tc1,vc1,tc2,vc2)


imprimeVetor(solexata)
print 'Soluçao Exata'

erroabsoluto = erroabsoluto(vetorMDF,solexata,npontos)

imprimeVetor(erroabsoluto)
print 'Erro absoluto'

plt.plot(abscissas,solnumerica,'bo-')   #gráfico soluçao numerica
    
plt.plot(abscissas,solexata,'xr-')   #gráfico solucao exata

#Identificaçao dos eixos
axes = plt.gca()
axes.set_xlabel('x')
axes.set_ylabel('u(x)')


plt.show()    #apresenta os gráficos plotados
