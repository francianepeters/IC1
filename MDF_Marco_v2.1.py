#!/usr/bin/env python
# *-* coding: utf-8 *-*
#https://www.wolframalpha.com/input/?i=u%3D%28%28%28x^4%29%2F12%29-%28%28x^3%29%2F6%29%29%2B%28%281%2F12%29%29x
#http://interactivepython.org/NWUbZ/LMdYZ/LpOMZ/courselib/static/pythonds/Graphs/graphintro.html

#Recado 25.06.2015
#A fazer:
#a) Criterio de parada;
#b) Corrigir solucoes exatas;
#Feito:
#a) Tamanho do dominio funciona para solucao numerica;
#b) Legenda no grafico
############################

import math
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
def MDFcoef1d(L,npontos,tc1,vc1,tc2,vc2):  
    h=(1.0*L)/(npontos-1) #h-dimensao dos trechos particionados

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
def vetorpoisson(k2,L,npontos,tc1,vc1,tc2,vc2):
    vetor=[0]*npontos
    h=(1.0*L)/(npontos-1)
    vetor[0]=vc1

    if k2 == 1:
        for i in range (1,npontos-1):
            x=i*h
            funcao=(x**2)-x
            vetor[i]=funcao
    if k2 == 2:
        for i in range (1,npontos-1):
            x=i*h
            funcao=x
            vetor[i]=funcao
    if k2 == 3:
        for i in range (1,npontos-1):
            x=i*h
            funcao=x*(x-1)*(x-0.5)
            vetor[i]=funcao
    if k2 == 4:
        for i in range (1,npontos-1):
            x=i*h
            funcao=math.exp(x)
            vetor[i]=funcao
    if k2 == 5:
        for i in range (1,npontos-1):
            x=i*h
            funcao=math.sin(x) + math.cos(x)
            vetor[i]=funcao
        
    vetor[npontos-1]=vc2

    return vetor

#Montar o vetor resposta U(x) para os pontos x desejados
#Laplace solucionado via analitica
def exatalaplace1d(L,npontos,tc1,vc1,tc2,vc2):  
    h=(1.0*L)/(npontos-1)
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
def exatapoisson1d(k2,L,npontos,tc1,vc1,tc2,vc2):  
    h=(1.0*L)/(npontos-1)
    exata=[0]*npontos
       
    if tc1==0 and tc2==0:
        if k2==1:
            for i in range (npontos):
                x=i*h
                u=(((x**4)/12.0)-((x**3)/6.0))+((1.0/12.0)+vc2-vc1)*x + vc1
                exata[i]=u
        if k2==2:
            for i in range (npontos):
                x=i*h
                u=((x**3)/6.0)-((1.0/6.0) + vc1 - vc2)*x + vc1
                exata[i]=u
        if k2==3:
            for i in range (npontos):
                x=i*h
                u=(((x**5)/20)-((x**4)/8)+((x**3)/12)) -(vc1-vc2+(1.0/120.0))*x + vc1
                exata[i]=u
        if k2==4:
            for i in range (npontos):
                x=i*h
                u=(math.exp(x)) - (math.e - 1 + vc1 - vc2)*x - 1+vc1
                exata[i]=u
        if k2==5:
            for i in range (npontos):
                x=i*h
                u=-(math.sin(x)+math.cos(x)) - (1+vc1-vc2-math.sin(1)-math.cos(1))*x +1+vc1
                exata[i]=u
            
    if tc1==0 and tc2==1:
        if k2==1:
            for i in range (npontos):
                x=i*h
                u=(((x**4)/12.0)-((x**3)/6.0))+(vc2 + (1/6.0))*x + vc1
                exata[i]=u
        if k2==2:
            for i in range (npontos):
                x=i*h
                u=((x**3)/6.0)-(0.5-vc2)*x + vc1
                exata[i]=u
        if k2==3:
            for i in range (npontos):
                x=i*h
                u=(((x**5)/20)-((x**4)/8)+((x**3)/12)) -((1.0/120.0)-vc2)*x + vc1
                exata[i]=u
        if k2==4:
            for i in range (npontos):
                x=i*h
                u=(math.exp(x)) - (math.e-vc2)*x - 1+vc1
                exata[i]=u
        if k2==5:
            for i in range (npontos):
                x=i*h
                u=-(math.sin(x)+math.cos(x)) - (math.sin(1)-math.cos(1)-vc2)*x +1+vc1
                exata[i]=u

    if tc1==1 and tc2==0:
        if k2==1:
            for i in range (npontos):
                x=i*h
                u=(((x**4)/12.0)-((x**3)/6.0))+(vc1)*x + ((1/12.0)+vc2-vc1)
                exata[i]=u
        if k2==2:
            for i in range (npontos):
                x=i*h
                u=((x**3)/6.0)+(vc1)*x -((1.0/6.0) + vc1 - vc2)
                exata[i]=u
        if k2==3:
            for i in range (npontos):
                x=i*h
                u=(((x**5)/20)-((x**4)/8)+((x**3)/12)) + (vc1)*x - (vc1-vc2+(1.0/120.0))
                exata[i]=u
        if k2==4:
            for i in range (npontos):
                x=i*h
                u=(math.exp(x)) - (1-vc1)*x - (math.e - 1 + vc1 - vc2)
                exata[i]=u
        if k2==5:
            for i in range (npontos):
                x=i*h
                u=-(math.sin(x)+math.cos(x)) - (1+vc1)*x - (1+vc1-vc2-math.sin(1)-math.cos(1))
                exata[i]=u

    return exata

def erroabsoluto(solnumerica,solexata,npontos):
    erro=[0]*npontos
    for i in range (npontos):
        erro[i]=solnumerica[i]-solexata[i]
    return erro

##############################################################################
#Estrutura principal do programa:
##############################################################################

k1=input('Digite "1" para resolver Laplace; "2" para Poisson:')

if k1==2:
    k2=input('Escolha uma das funçoes da lista, digitando seu respectivo numero:\n'
             '(1) => f(x) = x(x-1)\n'
             '(2) => f(x) = x\n'
             '(3) => f(x) = x(x-1)(x-0,5)\n'
             '(4) => f(x) = e**x\n'
             '(5) => f(x) = senx + cosx\n')
    
L=input('tamanho do dominio do problema:')
    
npontos=input('número de pontos para discretizacao:')

verificanpts(npontos)

tc1=input('tipo de condicao de contorno da esquerda:')
vc1=input('valor de condicao de contorno da esquerda:')

tc2=input('tipo de condicao de contorno da direita:')
vc2=input('valor de condicao de contorno da direita:')

verificacc(tc1,tc2)

matrizcoef=MDFcoef1d(L,npontos,tc1,vc1,tc2,vc2)

imprimeMatriz(matrizcoef)
print 'Matriz de Coeficientes MDF'

if k1==1:
    vetorMDF=vetorlaplace(npontos,tc1,vc1,tc2,vc2)

elif k1==2:
    vetorMDF=vetorpoisson(k2,L,npontos,tc1,vc1,tc2,vc2)

imprimeVetor(vetorMDF)
print 'Vetor MDF'

solnumerica=linalg.solve(matrizcoef,vetorMDF)

imprimeVetor(solnumerica)
print 'Soluçao Numerica'

abscissas=numpy.linspace(0,L,npontos)

imprimeVetor(abscissas)
print 'Pontos discretizados'

if k1==1:
    solexata=exatalaplace1d(L,npontos,tc1,vc1,tc2,vc2)

elif k1==2:
    solexata=exatapoisson1d(k2,L,npontos,tc1,vc1,tc2,vc2)


imprimeVetor(solexata)
print 'Soluçao Exata'

erroabsoluto = erroabsoluto(solnumerica,solexata,npontos)

imprimeVetor(erroabsoluto)
print 'Erro absoluto'

plt.plot(abscissas,solnumerica,'bo-', label= 'Numerica')   #gráfico soluçao numerica
    
plt.plot(abscissas,solexata,'xr-', label='Exata')   #gráfico solucao exata

plt.legend(loc='upper left')

#Identificaçao dos eixos
axes = plt.gca()
axes.set_xlabel('x')
axes.set_ylabel('u(x)')


plt.show()    #apresenta os gráficos plotados
