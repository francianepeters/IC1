#!/usr/bin/env python
# *-* coding: utf-8 *-*
import math
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
def montavetor(npontos,tc1,vc1,tc2,vc2,problema):
    vetor=[0]*npontos

    vetor[0]=vc1
    vetor[npontos-1]=vc2

    h=1.0/(npontos-1)

    print problema
    print h

#Mudança: antes o "for" era (1,npontos-1)
    if problema==1:
        for i in range (1,npontos-2):
            x=i*h
            vetor[i]=x
        

    return vetor
#Problema de Poisson para f(x) = x² - x
def vetorpoisson(k2,npontos,tc1,vc1,tc2,vc2):
    vetor=[0]*npontos
    h=1.0/(npontos-1)
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
#Poisson solucionado via analitica para f(x) = x² - x
def exatapoisson1d(k2,npontos,tc1,vc1,tc2,vc2):  
    h=1.0/(npontos-1)
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
        if k2==6:
            print "Nao foi calculado"
                
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
        if k2==6:
            print "Nao foi calculado"       

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
        if k2==6:
            x=i*h
            u=((3*x**10)/90)+((4*x**8)/56)+((45*x**2)/2)+vc1+vc2-((3*x**10)/90)-((4*x**8)/56)-((45*x**2)/2)-vc1x
            exalta[i]=u
            
    return exata

def erroabsoluto(solnumerica,solexata,npontos):
    erro=[0]*npontos
    for i in range (npontos):
        erro[i]=solnumerica[i]-solexata[i]
    return erro


##############################################################################
#Programa principal 
##############################################################################
problema=input('digite 1 para laplace e 2 para poisson:')
if problema==2:
    k2=input('Escolha uma das funçoes da lista, digitando seu respectivo numero:\n'
             '(1) => f(x) = x(x-1)\n'
             '(2) => f(x) = x\n'
             '(3) => f(x) = x(x-1)(x-0,5)\n'
             '(4) => f(x) = e**x\n'
             '(5) => f(x) = senx + cosx\n')
             '(6) => f(x) = ((3x**8)+(4x**6)+45')
    

npontos=input('número de pontos para discretizacao:')

#if(npontos<=1)then
#('número de pontos incorreto')

tc1=input('tipo de condicao de contorno da esquerda:')
vc1=input('valor de condicao de contorno da esquerda:')


tc2=input('tipo de condicao de contorno da direita:')
vc2=input('valor de condicao de contorno da direita:')


matriz=laplace1d(npontos,tc1,vc1,tc2,vc2)

imprimiMatriz(matriz)
print 'Matriz de Coeficientes MDF'

if problema==1:
    vetor=vetorlaplace(npontos,tc1,vc1,tc2,vc2)

elif problema==2:
    vetor=vetorpoisson(k2,npontos,tc1,vc1,tc2,vc2)

imprimiVetor(vetor)

solucao=linalg.solve(matriz,vetor)

imprimiVetor(solucao)

x=numpy.linspace(0,1,npontos)

imprimiVetor(x)


if problema==1:
   exata=exatalaplace1d(npontos,tc1,vc1,tc2,vc2)
elif problema==2:
    exata=exatapoisson1d(k2,npontos,tc1,vc1,tc2,vc2)

imprimiVetor(exata)


plt.plot(x,solucao,'bo-')   # estou a fazer o gráfico dos dois vectores.

    
plt.plot(x,exata,'rx-')   # estou a fazer o gráfico dos dois vectores.


plt.show()    # Esta linha de código serve para mostrar o gráfico total.
